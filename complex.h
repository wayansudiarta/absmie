#ifndef _NR_COMPLEX_H_
#define _NR_COMPLEX_H_

#ifndef _FCOMPLEX_DECLARE_T_
typedef struct FCOMPLEX {float r,i;} fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

fcomplex Cadd(fcomplex a, fcomplex b);
fcomplex Csub(fcomplex a, fcomplex b);
fcomplex Cmul(fcomplex a, fcomplex b);
fcomplex Complex(float re, float im);
fcomplex Conjg(fcomplex z);
fcomplex Cdiv(fcomplex a, fcomplex b);
float Cabs(fcomplex z);
fcomplex Csqrt(fcomplex z);
fcomplex RCmul(float x, fcomplex a);

/* this is for complex sine and cosine */
fcomplex Csin(fcomplex z);
fcomplex Ccos(fcomplex z);

#else /* ANSI */
/* traditional - K&R */

fcomplex Cadd();
fcomplex Csub();
fcomplex Cmul();
fcomplex Complex();
fcomplex Conjg();
fcomplex Cdiv();
float Cabs();
fcomplex Csqrt();
fcomplex RCmul();
fcomplex Ccos();
fcomplex Csin();

#endif /* ANSI */

#endif /* _NR_COMPLEX_H_ */

