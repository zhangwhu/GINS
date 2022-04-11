#ifndef PPPSTATE_H
#define PPPSTATE_H

/* number and index of states */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)				//频率
#define NP(opt)     ((opt)->dynamics?9:3)									//位置
#define NC(opt)     (NSYS)														//系统
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))//对流层
#define NI(opt)     ((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)					//电离层（每颗卫星一量）
#define NS(opt)		((opt)->nf>2?MAXSAT:0)									//IFB
#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt)+NS(opt))				//P+C+trop+Iono+IFB
#define NB(opt)     (NF(opt)*MAXSAT)										//模糊度
#define NX(opt)     (NR(opt)+NB(opt))										//总状态
#define IC(s,opt)   (NP(opt)+(s))
#define IT(opt)     (NP(opt)+NC(opt))
#define II(s,opt)   (NP(opt)+NC(opt)+NT(opt)+(s)-1)
#define IS(s,opt)	(NP(opt)+NC(opt)+NT(opt)+NI(opt)+(s)-1)
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)

#define REL_HUMI    0.5         /* relative humidity for saastamoinen model */

#endif
