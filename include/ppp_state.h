#ifndef PPPSTATE_H
#define PPPSTATE_H

/* number and index of states */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)				
#define NP(opt)     (9)	                                                    //状态有必要设置为9维吗？								
#define NC(opt)     (NSYS)													
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))
#define NI(opt)     ((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)					//����㣨ÿ������һ����
#define NS(opt)		((opt)->nf>2?MAXSAT:0)									//IFB
#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt)+NS(opt))				//P+C+trop+Iono+IFB
#define NB(opt)     (NF(opt)*MAXSAT)										//ģ����
#define NX(opt)     (NR(opt)+NB(opt))										//��״̬
#define IC(s,opt)   (NP(opt)+(s))
#define IT(opt)     (NP(opt)+NC(opt))
#define II(s,opt)   (NP(opt)+NC(opt)+NT(opt)+(s)-1)
#define IS(s,opt)	(NP(opt)+NC(opt)+NT(opt)+NI(opt)+(s)-1)
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)

#define REL_HUMI    0.5         /* relative humidity for saastamoinen model */

#endif
