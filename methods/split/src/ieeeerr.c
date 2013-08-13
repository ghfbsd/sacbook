/* Function to be error handler for IEEE errors.  The handler arranges that
   an ABRT signal is generated at the exception point so that debuggers may 
   be used to figure out where and why the exception occurred. 

   Usage:  Callable from FORTRAN as
      call ieeeset(what)

      what - character parameter which may be specifically named
	 exception types (taken from SunOS 4.1 documentation of
	 'ieee_handler(3)'), or which may be the value 'environment'.
	 This describes the type of exceptions that will trigger an
	 error message and abort the running program.
	 
	 If the 'environment' is given, the value used is the value
	 that the environment variable IEEE_TRAP has.  This allows
	 trapping to be enabled on a run-by-run basis.

	 Types:
	    "inexact" - inexact result
	    "division" - division by zero
	    "underflow" - result of operation too small to represent
	    "overflow" - result of operation too large to represent
	    "invalid" - invalid operand
	    "all" - synonym for all of above
	    "common" - synonym for division, overflow and inexact

   G. Helffrich CIW/DTM Jul. 25 1991

   Modified for Solaris/SysV use - G. Helffrich Autumn 1997

   Modified for DEC Alpha use - G. Helffrich 19 Oct. 1998

   Modified for Linux/i386 use - G. Helffrich 18 June 1999

   Modified for Linux/sparc use - G. Helffrich 5 June 2000

   Modified for MacOS X/ppc use - G. Helffrich 3 June 2003, 4 Aug. 2003

   Modified for FreeBSD/i386 use - G. Helffrich 4 Aug. 2003
 
--------------------------------------------------------------------------------
  Copyright (c) 2013 by G. Helffrich.
  All rights reserved.
 
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution.  
    * The names of the authors may not be used to endorse or promote products
      derived from this software without specific prior written permission.
  THERE IS NO WARRANTY EXPRESSED OR IMPLIED, NOR LIABILITY FOR DAMAGES ASSUMED,
  IN PROVIDING THIS SOFTWARE.
--------------------------------------------------------------------------------
*/

#include <signal.h>
#include <stdio.h>
#include <strings.h>
#if defined (_ALPHA_SIGNAL_H_)
#include <machine/fpu.h>
#elif defined (__ppc__)
#include <architecture/ppc/fp_regs.h>
#include <mach/mach.h>
#include <pthread.h>
#define _FPU_MASK_IM (1<<7)
#define _FPU_MASK_OM (1<<6)
#define _FPU_MASK_UM (1<<5)
#define _FPU_MASK_ZM (1<<4)
#define _FPU_MASK_PM (1<<3)
#define sc_pc sc_ir
#elif defined (__linux__)
/*  Linux variant of FPU control.  Danger, icebergs ahead.  */
#include <fpu_control.h>
#if defined(_SPARC_FPU_CONTROL_H)
/*  Linux 2.0.35 kernel code doesn't define these, we have to ourselves.  */

#define CONTEXT_PC context->sigc_pc
#define CONTEXT_FP_SWORD 0 /*dummy, is context->fpstate->sw for I386 */
#define CONTEXT_FP_CWORD 0 /*dummy, context->fpstate->cw for I386 */

/*  2.0.35 kernel headers don't provide access to FPU status word in the
    signal context, so these defines turn ieeeset into a no-op on SPARC
    G. Helffrich/5 June 2000 */

#ifndef _FPU_MASK_IM
#define _FPU_MASK_IM 0x00
#endif
#ifndef _FPU_MASK_DM
#define _FPU_MASK_DM 0x00
#endif
#ifndef _FPU_MASK_ZM
#define _FPU_MASK_ZM 0x00
#endif
#ifndef _FPU_MASK_OM
#define _FPU_MASK_OM 0x00
#endif
#ifndef _FPU_MASK_UM
#define _FPU_MASK_UM 0x00
#endif
#ifndef _FPU_MASK_PM
#define _FPU_MASK_PM 0x00
#endif
#endif
#if defined(_I386_FPU_CONTROL_H)
#define CONTEXT_PC context->eip
#define CONTEXT_FP_SWORD context->fpstate->sw
#define CONTEXT_FP_CWORD context->fpstate->cw
#endif
#else
#include <floatingpoint.h>
/* SunOS 5.7 compilers without ucbcc need help identifying themselves */
#ifdef _SYS_IEEEFP_H
#define SYSV
#undef BSD
#endif
#endif
#ifndef NULL
#define NULL 0
#endif
/*
 * This trick is used to distinguish between SYSV and V7 systems.
 * We assume that L_ctermid is only defined in stdio.h in SYSV
 * systems, but not in V7 or Berkeley UNIX.
 */
#if defined(L_ctermid) && !defined(_SYS_IEEEFP_H)
#ifdef SYSV
#undef SYSV
#endif
#define BSD
#endif
/* Some V7 systems define       L_ctermid - we list those here */
#if defined(BSD) && !defined(_SYS_IEEEFP_H)
#undef SYSV
#endif
#ifdef __linux__
#undef BSD
#undef SYSV
#endif

struct ieee_gaz {
   int ieee_code;
   char *ieee_words;
} ieee_decode[] = {
/* BSD/SunOS 4.x */
#ifdef FPE_INTOVF_TRAP
   {FPE_INTOVF_TRAP,	"integer overflow"},
#endif
#ifdef FPE_INTDIV_TRAP
   {FPE_INTDIV_TRAP,	"integer divide by zero"},
#endif
#ifdef FPE_FLTINEX_TRAP
   {FPE_FLTINEX_TRAP,	"floating inexact result"},
#endif
#ifdef FPE_FLTDIV_TRAP
   {FPE_FLTDIV_TRAP,	"floating divide by zero"},
#endif
#ifdef FPE_FLTUND_TRAP
   {FPE_FLTUND_TRAP,	"floating underflow"},
#endif
#ifdef FPE_FLTOPERR_TRAP
   {FPE_FLTOPERR_TRAP,	"floating operand error"},
#endif
#ifdef FPE_FLTOVF_TRAP
   {FPE_FLTOVF_TRAP,	"floating overflow"},
#endif
#ifdef FPE_CHKINST_TRAP
   {FPE_CHKINST_TRAP,   "CHK [CHK2] instruction"},
#endif
#ifdef FPE_TRAPV_TRAP
   {FPE_TRAPV_TRAP,     "TRAPV [cpTRAPcc TRAPcc] instr"},
#endif
#ifdef FPE_FLTBSUN_TRAP
   {FPE_FLTBSUN_TRAP,   "branch or set on unordered cond"},
#endif
#ifdef FPE_FLTNAN_TRAP
   {FPE_FLTNAN_TRAP,    "floating Not-A-Number"},
#endif
/* SYSV/SunOS 5.x */
#ifdef FPE_FLTINV
   {FPE_FLTINV, "invalid operand"},
#endif
#ifdef FPE_FLTRES
   {FPE_FLTRES, "inexact"},
#endif
#ifdef FPE_FLTDIV
   {FPE_FLTDIV, "division-by-zero"},
#endif
#ifdef FPE_FLTUND
   {FPE_FLTUND, "underflow"},
#endif
#ifdef FPE_FLTOVF
   {FPE_FLTUND, "overflow"},
#endif
/* Digital UNIX on DEC Alpha */
#ifdef IEEE_STATUS_INV
   {IEEE_STATUS_INV, "Invalid operation"},
#endif
#ifdef IEEE_STATUS_DZE
   {IEEE_STATUS_DZE, "Divide by 0"},
#endif
#ifdef IEEE_STATUS_OVF
   {IEEE_STATUS_OVF, "Overflow"},
#endif
#ifdef IEEE_STATUS_UNF
   {IEEE_STATUS_UNF, "Underflow"},
#endif
#ifdef IEEE_STATUS_INE
   {IEEE_STATUS_INE, "Inexact"},
#endif
/* Linux i386/fpu_control.h */
#ifdef _FPU_MASK_IM
   {_FPU_MASK_IM, "Invalid"},
#endif
#ifdef _FPU_MASK_DM
   {_FPU_MASK_DM, "Denormalized"},
#endif
#ifdef _FPU_MASK_ZM
   {_FPU_MASK_ZM, "Divide by 0"},
#endif
#ifdef _FPU_MASK_OM
   {_FPU_MASK_OM, "Overflow"},
#endif
#ifdef _FPU_MASK_UM
   {_FPU_MASK_UM, "Underflow"},
#endif
#ifdef _FPU_MASK_PM
   {_FPU_MASK_PM, "Inexact"},
#endif
#if defined(__linux__) & defined(__i386__)
#define _LINUXI386_OVR 0x10000
   {_LINUXI386_OVR, "coprocessor segment overrun"},
#define _LINUXI386_OVF 0x00040
   {_LINUXI386_OVR, "coprocessor stack overflow"},
#endif
   { -1, NULL}
};

#if defined(__MACH__)
/* These functions change the machine registers but Mach kernel will change
   them back based on its internal state, so they aren't very useful.  Here
   for reference, however.
*/
#if 0

void setfpucw(unsigned int v)
{
   union {
     unsigned int d[2];
     double align;
   } u;
   u.d[0] = 0; u.d[1] = v;
   asm ("lfd f0,%0\n\tmtfsf 0xff,f0" : : "m" (u.d[0]) : "f0" );
}

unsigned int getfpucw(void)
{
   union {
     unsigned int d[2];
     double align;
   } u;
   asm ("mffs f0\n\tstfd f0,%0" : "=m" (u.d[0]) : : "f0" );
   return u.d[1];
}
#endif

/* Many thanks to the yorick-1.5 developers who published this code - I couldn't
   have figured this out without their help.

   G. Helffrich 4 Aug. 2003
*/

#define FE0_MASK (1<<11)
#define FE1_MASK (1<<8)
/* FE0  FE1   exceptions enabled if either FE0 or FE1 set
 *  0    0    -- floating-point exceptions disabled
 *  0    1    -- floating-point imprecise nonrecoverable
 *  1    0    -- floating-point imprecise recoverable
 *  1    1    -- floating-point precise mode
 */

/* a thread cannot get or set its own MSR bits */
static void *
fpu_fpe_enable(void *arg)  
{
  thread_t t = *(thread_t *)arg;
  struct ppc_thread_state state;
  unsigned int state_size = PPC_THREAD_STATE_COUNT;
  if (thread_get_state(t, PPC_THREAD_STATE,
                       (natural_t *)&state, &state_size) == KERN_SUCCESS) {
    state.srr1 |= FE1_MASK;      
    thread_set_state(t, PPC_THREAD_STATE, (natural_t *)&state, state_size);
  }
  return 0;
}
#endif

void
ieeeerr(code,where)
   int code;
   char *where;
{
   struct ieee_gaz *sp;
   char *words = "(unknown code)";

   /* Turn code into words, then print message. */
   for (sp=ieee_decode; sp->ieee_code != -1; sp++)
      if (sp->ieee_code == code) {
	 words = sp->ieee_words; break;
      }
   fprintf(stderr,"IEEE error %x at %X: %s.\n",
      code,where,words);

   /* Block abort signal, send an abort signal to ourselves,
      and rely on restored signal mask on return from signal
      handler to unblock the abort signal and stop the program. */
#ifdef SYSV
   (void) sighold(SIGABRT);
#else
   (void) sigblock(sigmask(SIGABRT));
#endif
   (void) kill(getpid(),SIGABRT);
}

#ifdef SYSV
#include <siginfo.h>
#include <ucontext.h>

void
ieeeexf(sig, sip, scp)
   int sig;
   siginfo_t *sip;
   ucontext_t *scp;
{
   fpregset_t *uc = &scp->uc_mcontext.fpregs;
   int code = sip->si_code;

   ieeeerr(code,scp->uc_mcontext.gregs[REG_PC]);
}
#endif

#if defined(BSD) & !defined(_ALPHA_SIGNAL_H_)

void
ieeeexf(sig, code, scp, addr)
   int sig, code;
   struct sigcontext *scp;
   char *addr;
{
   ieeeerr(code,scp->sc_pc);
}
#endif

#if defined(BSD) & defined(_ALPHA_SIGNAL_H_)
#include <ucontext.h>
/* The pc isn't set right by this, must be incorrect returned from getcontext.
   Any ideas on how to get the right pc on a DEC Alpha would be appreciated.
   G. Helffrich 19 June 1999
*/

void
ieeeexf(sig)
   int sig;
{
   unsigned long cntl, fpcr, pc;
   ucontext_t context;

   getcontext(&context);
   cntl = context.uc_mcontext.sc_fpcr; pc = context.uc_mcontext.sc_pc;
   ieeeerr((cntl & IEEE_STATUS_MASK) >> IEEE_STATUS_TO_EXCSUM_SHIFT,pc);
}
#endif

#ifdef __linux__
#include <asm/sigcontext.h>
/* fexcp package by W. Metzenthen helped guide development of Linux/386
   section.

   Should be a better way to get the prevailing FP status other than to
   send us an interrupt to get a peek at the context.

   G. Helffrich 19 June 1999
*/

static int fpstat;

static void
getfpexf(sig)
   int sig;
{  
   /* Resort to low cunning to get pointer to context following last
      parameter in stack */
   int code;
   char *addr;
   struct sigcontext_struct *context = (void *)(&sig + 1);

   if (sig != SIGUSR1)
      fprintf(stderr,
         "**IEEERR:  Unexpected signal %d, wanted %d (SIGUSR1)\n",sig,SIGUSR1);
   fpstat = CONTEXT_FP_CWORD & 0xffff;
}

unsigned short
getfpstatus() {
	signal(SIGUSR1,getfpexf);
	raise(SIGUSR1);
	return ((unsigned short)fpstat);
}

static void
ieeeexf(sig)
   int sig;
{  
   /* Resort to low cunning to get pointer to context following last
      parameter in stack */
   int code;
   char *addr;
   struct sigcontext_struct *context = (void *)(&sig + 1);

   if (sig != SIGFPE) {
      fprintf(stderr,
         "**IEEERR:  Unexpected signal %d, wanted %d (SIGFPE)!\n",sig,SIGFPE);
      signal(SIGFPE,ieeeexf);
   } else switch (context->trapno) {
      case  0: /* Normal, integer division error */
      case 16: /* Normal, floating point error */
	 code = CONTEXT_FP_SWORD & 0x007f;
         /* Clear imprecise bit if set in conjunction with anything else */
	 if ((code & ~_FPU_MASK_PM) && (code & _FPU_MASK_PM))
	    code &= ~_FPU_MASK_PM;
         /* Report stack overflow if set in conjunction with invalid op */
#if defined(_LINUXI386_OVF)
	 if ((code & _FPU_MASK_IM) && (code & _LINUXI386_OVF))
	    code = _LINUXI386_OVF;
#endif
	 ieeeerr(code, CONTEXT_PC);
	 break;
#if defined(__i386__)
      case 9: /* Coprocessor segment overrun */
         ieeeerr(_LINUXI386_OVR, context->eip);
	 break;
#endif
      default:
	 fprintf(stderr,
            "**IEEERR:  SIGFPE trap %d ignored.\n",context->trapno);
	 signal(SIGFPE,ieeeexf);
      }
}
#endif

ieeeset_(what)
   char *what;
{
   char *ptr, *getenv();

   if (0 == strncmp(what,"environment",11)) {
      /* Get default value from environment */
      if (NULL == (ptr = getenv("IEEE_TRAP"))) return;
   } else
      ptr = what;
      
#if defined (_ALPHA_SIGNAL_H_)
   { unsigned long mask = 0;
     if (strstr(ptr,"inexact")) mask |= IEEE_TRAP_ENABLE_INE;
     if (strstr(ptr,"division")) mask |= IEEE_TRAP_ENABLE_DZE;
     if (strstr(ptr,"underflow")) mask |= IEEE_TRAP_ENABLE_UNF;
     if (strstr(ptr,"overflow")) mask |= IEEE_TRAP_ENABLE_OVF;
     if (strstr(ptr,"invalid")) mask |= IEEE_TRAP_ENABLE_INV;
     if (strstr(ptr,"all")) mask |= IEEE_TRAP_ENABLE_INE
				 | IEEE_TRAP_ENABLE_DZE
				 | IEEE_TRAP_ENABLE_UNF
				 | IEEE_TRAP_ENABLE_OVF
				 | IEEE_TRAP_ENABLE_INV;
     if (strstr(ptr,"common")) mask |= IEEE_TRAP_ENABLE_DZE
				 | IEEE_TRAP_ENABLE_INV
				 | IEEE_TRAP_ENABLE_OVF;
     if (mask) {
	ieee_set_fp_control(mask);
	/* Call signal handler */
	signal(SIGFPE,ieeeexf);
     }
   }
#elif defined (__linux__)
   {
     unsigned short mask = 0;
     unsigned short ctrl;

     if (strstr(ptr,"inexact")) mask |= _FPU_MASK_PM ;
     if (strstr(ptr,"division")) mask |= _FPU_MASK_ZM ;
     if (strstr(ptr,"underflow")) mask |= _FPU_MASK_UM ;
     if (strstr(ptr,"overflow")) mask |= _FPU_MASK_OM ;
     if (strstr(ptr,"invalid")) mask |= _FPU_MASK_IM ;
     if (strstr(ptr,"all")) mask |= _FPU_MASK_PM 
				 | _FPU_MASK_ZM
                                 | _FPU_MASK_UM
                                 | _FPU_MASK_OM
                                 | _FPU_MASK_IM ;
     if (strstr(ptr,"common")) mask |= _FPU_MASK_ZM
                                    | _FPU_MASK_OM
                                    | _FPU_MASK_IM ;
     ctrl = getfpstatus();
     if (!ctrl) __setfpucw(ctrl=_FPU_IEEE);
     if (mask) {
        /* Code is Linux/i386 specific */
        ctrl &= ~mask; __setfpucw (ctrl);
	/* Call signal handler */
	signal(SIGFPE,ieeeexf);
     }
   }
#elif defined (__ppc__)
   {
     unsigned mask = 0;
     unsigned ctrl;
     pthread_t enabler;
     thread_t self = mach_thread_self();

     if (strstr(ptr,"inexact")) mask |= _FPU_MASK_PM ;
     if (strstr(ptr,"division")) mask |= _FPU_MASK_ZM ;
     if (strstr(ptr,"underflow")) mask |= _FPU_MASK_UM ;
     if (strstr(ptr,"overflow")) mask |= _FPU_MASK_OM ;
     if (strstr(ptr,"invalid")) mask |= _FPU_MASK_IM ;
     if (strstr(ptr,"all")) mask |= _FPU_MASK_PM 
				 | _FPU_MASK_ZM
                                 | _FPU_MASK_UM
                                 | _FPU_MASK_OM
                                 | _FPU_MASK_IM ;
     if (strstr(ptr,"common")) mask |= _FPU_MASK_ZM
                                    | _FPU_MASK_OM
                                    | _FPU_MASK_IM ;
     if (mask) {
#if 0
        /* Clear exception enable masks */
	ctrl = getfpucw () & ~( _FPU_MASK_PM | _FPU_MASK_ZM | _FPU_MASK_UM |
                               _FPU_MASK_OM | _FPU_MASK_IM );
        setfpucw (ctrl & mask);

#else
        ppc_fp_scr_t r = get_fp_scr();     

	/* turn off exception bits to start with none set */
	r.fx = r.fex = r.ox = r.ux = r.zx = r.xx =
	   r.vx_zdz = r.vx_imz = r.vx_xvc = r.vx_cvi = r.vx_soft = 0;

	/* Set desired exception bits */
	if (mask & _FPU_MASK_PM) r.xe = 1;
	if (mask & _FPU_MASK_ZM) r.ze = 1;
	if (mask & _FPU_MASK_UM) r.ue = 1;
	if (mask & _FPU_MASK_OM) r.oe = 1;
	if (mask & _FPU_MASK_IM) r.ve = 1;
	set_fp_scr(r);
#endif

	/* Enable FP exceptions in our thread */
	if (!pthread_create(&enabler, 0, fpu_fpe_enable, &self))
	   pthread_join(enabler, 0);

	/* Call signal handler */
  	signal(SIGFPE,ieeeexf);
     }
   }
#elif defined (BSD)
   {
     fp_except_t mask = 0;
     fp_except_t ctrl;

     if (strstr(ptr,"inexact")) mask |= FP_X_IMP ;
     if (strstr(ptr,"division")) mask |= FP_X_DZ ;
     if (strstr(ptr,"underflow")) mask |= FP_X_UFL ;
     if (strstr(ptr,"overflow")) mask |= FP_X_OFL ;
     if (strstr(ptr,"invalid")) mask |= FP_X_INV ;
     if (strstr(ptr,"all")) mask |= FP_X_IMP 
				 | FP_X_DZ
                                 | FP_X_UFL
                                 | FP_X_OFL
                                 | FP_X_INV ;
     if (strstr(ptr,"common")) mask |= FP_X_DZ
                                    | FP_X_OFL
                                    | FP_X_INV ;
     if (mask) {
        /* Clear exception enable masks */
	ctrl = fpgetmask() & ~( FP_X_IMP | FP_X_DZ | FP_X_UFL | FP_X_OFL
                                 | FP_X_INV ) ;
        (void) fpsetmask (ctrl | mask);
	/* Call signal handler */
  	signal(SIGFPE,ieeeexf);
     }
   }
#else
   if (0 != ieee_handler("set",ptr,ieeeexf))
      fprintf(stderr,"IEEESET: Error calling handler with '%s'\n",
      ptr);
#endif
}
