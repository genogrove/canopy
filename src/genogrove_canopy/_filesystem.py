# SPDX-License-Identifier: GPL-3.0-or-later
"""OS filesystem policy installed in the child after trusted imports, before query code."""


def restrict_filesystem(roots):
    """Self-contained for source injection. Never fall back to Python-only isolation."""
    import ctypes
    import json
    import os
    import sys

    if sys.platform == "darwin":
        lib = ctypes.CDLL("/usr/lib/libsandbox.dylib", use_errno=True)
        lib.sandbox_init.argtypes = [ctypes.c_char_p, ctypes.c_uint64,
                                     ctypes.POINTER(ctypes.c_char_p)]
        lib.sandbox_init.restype = ctypes.c_int
        lib.sandbox_free_error.argtypes = [ctypes.c_char_p]
        # Apply after imports so neither the interpreter nor its native libraries need
        # broad read grants. Metadata lookups remain allowed; file contents do not.
        profile = '(version 1)(allow default)(deny file-read-data)(deny file-write*)'
        for path in roots:
            kind = "subpath" if os.path.isdir(path) else "literal"
            profile += '(allow file-read-data (%s %s))' % (kind, json.dumps(path, ensure_ascii=False))
        error = ctypes.c_char_p()
        if lib.sandbox_init(profile.encode(), 0, ctypes.byref(error)):
            message = error.value.decode() if error.value else "sandbox_init failed"
            lib.sandbox_free_error(error)
            raise RuntimeError("Filesystem isolation unavailable: " + message)
    elif sys.platform == "linux" and os.uname().machine in ("x86_64", "aarch64"):
        # Landlock syscall numbers on these architectures; ABI 3 is required to deny
        # O_RDONLY|O_TRUNC as well as write opens. See https://www.kernel.org/doc/html/latest/userspace-api/landlock.html.
        lib = ctypes.CDLL(None, use_errno=True)
        lib.syscall.restype = ctypes.c_long

        def checked(number, *args):
            result = lib.syscall(number, *args)
            if result < 0:
                raise OSError(ctypes.get_errno(), "Filesystem isolation unavailable (Landlock)")
            return result

        abi = checked(444, 0, 0, 1)  # landlock_create_ruleset: query ABI
        if abi < 3:
            raise RuntimeError("Filesystem isolation requires Landlock ABI 3 or newer")

        class Ruleset(ctypes.Structure):
            _fields_ = [("handled_access_fs", ctypes.c_uint64)]

        class PathBeneath(ctypes.Structure):
            _pack_ = 1
            _fields_ = [("allowed_access", ctypes.c_uint64), ("parent_fd", ctypes.c_int32)]

        rights = (1 << 15) - 1  # ABI 3: execute/read/write/create/remove/refer/truncate
        ruleset = Ruleset(rights)
        fd = checked(444, ctypes.byref(ruleset), ctypes.sizeof(ruleset), 0)
        try:
            for path in roots:
                try:
                    parent = os.open(path, os.O_PATH | os.O_CLOEXEC)
                except FileNotFoundError:
                    continue  # absent paths grant nothing; directories must exist at startup
                try:
                    read = (1 << 2) | ((1 << 3) if os.path.isdir(path) else 0)
                    rule = PathBeneath(read, parent)
                    checked(445, fd, 1, ctypes.byref(rule), 0)  # landlock_add_rule
                finally:
                    os.close(parent)
            if lib.prctl(38, 1, 0, 0, 0):  # PR_SET_NO_NEW_PRIVS
                raise OSError(ctypes.get_errno(), "Cannot set no_new_privs")
            checked(446, fd, 0)  # landlock_restrict_self: inherited and irreversible
        finally:
            os.close(fd)
    else:
        raise RuntimeError("Filesystem isolation requires macOS or Linux with Landlock ABI 3+")
