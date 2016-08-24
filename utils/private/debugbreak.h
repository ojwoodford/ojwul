/* Copyright (c) 2011-2015, Scott Tsai
 * 
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef DEBUG_BREAK_H
#define DEBUG_BREAK_H

#ifdef _MSC_VER

#include <Windows.h>
#define debug_break() if(IsDebuggerPresent()){__debugbreak();}

#else

#include <signal.h>
#include <stdbool.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/sysctl.h>

__inline__ static bool IsDebuggerPresent(void)
    // Returns true if the current process is being debugged (either 
    // running under the debugger or has a debugger attached post facto).
{
    int                 junk;
    int                 mib[4];
    struct kinfo_proc   info;
    size_t              size;

    // Initialize the flags so that, if sysctl fails for some bizarre 
    // reason, we get a predictable result.
    info.kp_proc.p_flag = 0;

    // Initialize mib, which tells sysctl the info we want, in this case
    // we're looking for information about a specific process ID.
    mib[0] = CTL_KERN;
    mib[1] = KERN_PROC;
    mib[2] = KERN_PROC_PID;
    mib[3] = getpid();

    // Call sysctl.
    size = sizeof(info);
    junk = sysctl(mib, sizeof(mib) / sizeof(*mib), &info, &size, NULL, 0);

    // We're being debugged if the P_TRACED flag is set.
    return ( (info.kp_proc.p_flag & P_TRACED) != 0 );
}

#ifdef __cplusplus
extern "C" {
#endif

enum {
	/* gcc optimizers consider code after __builtin_trap() dead.
	 * Making __builtin_trap() unsuitable for breaking into the debugger */
	DEBUG_BREAK_PREFER_BUILTIN_TRAP_TO_SIGTRAP = 0,
};

#if defined(__i386__) || defined(__x86_64__)
enum { HAVE_TRAP_INSTRUCTION = 1, };
__attribute__((gnu_inline, always_inline))
__inline__ static void trap_instruction(void)
{
	__asm__ volatile("int $0x03");
}
#elif defined(__thumb__)
enum { HAVE_TRAP_INSTRUCTION = 1, };
/* FIXME: handle __THUMB_INTERWORK__ */
__attribute__((gnu_inline, always_inline))
__inline__ static void trap_instruction(void)
{
	/* See 'arm-linux-tdep.c' in GDB source.
	 * Both instruction sequences below work. */
#if 1
	/* 'eabi_linux_thumb_le_breakpoint' */
	__asm__ volatile(".inst 0xde01");
#else
	/* 'eabi_linux_thumb2_le_breakpoint' */
	__asm__ volatile(".inst.w 0xf7f0a000");
#endif

	/* Known problem:
	 * After a breakpoint hit, can't stepi, step, or continue in GDB.
	 * 'step' stuck on the same instruction.
	 *
	 * Workaround: a new GDB command,
	 * 'debugbreak-step' is defined in debugbreak-gdb.py
	 * that does:
	 * (gdb) set $instruction_len = 2
	 * (gdb) tbreak *($pc + $instruction_len)
	 * (gdb) jump   *($pc + $instruction_len)
	 */
}
#elif defined(__arm__) && !defined(__thumb__)
enum { HAVE_TRAP_INSTRUCTION = 1, };
__attribute__((gnu_inline, always_inline))
__inline__ static void trap_instruction(void)
{
	/* See 'arm-linux-tdep.c' in GDB source,
	 * 'eabi_linux_arm_le_breakpoint' */
	__asm__ volatile(".inst 0xe7f001f0");
	/* Has same known problem and workaround
	 * as Thumb mode */
}
#elif defined(__aarch64__)
enum { HAVE_TRAP_INSTRUCTION = 1, };
__attribute__((gnu_inline, always_inline))
__inline__ static void trap_instruction(void)
{
	/* See 'aarch64-tdep.c' in GDB source,
	 * 'aarch64_default_breakpoint' */
	__asm__ volatile(".inst 0xd4200000");
}
#else
enum { HAVE_TRAP_INSTRUCTION = 0, };
#endif

__attribute__((gnu_inline, always_inline))
__inline__ static void debug_break(void)
{
    /* Only stop if debugger is attached. */
    if (!IsDebuggerPresent())
        return;
    
	if (HAVE_TRAP_INSTRUCTION) {
		trap_instruction();
	} else if (DEBUG_BREAK_PREFER_BUILTIN_TRAP_TO_SIGTRAP) {
		 /* raises SIGILL on Linux x86{,-64}, to continue in gdb:
		  * (gdb) handle SIGILL stop nopass
		  * */
		__builtin_trap();
	} else {
		#ifdef _WIN32
		/* SIGTRAP available only on POSIX-compliant operating systems
		 * use builtin trap instead */
		__builtin_trap();
		#else
		raise(SIGTRAP);
		#endif
	}
}

#ifdef __cplusplus
}
#endif

#endif

#endif