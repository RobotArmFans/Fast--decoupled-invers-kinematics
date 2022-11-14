#ifndef __INVER_KINEMATIC_H
#define __INVER_KINEMATIC_H 			   
#include "stm32f10x.h"

	
	typedef struct      //一个结构体
{
  float jion_1;
  float jion_2;
	float jion_3;
	float ola_Z1;
  float ola_Y;
	float ola_Z2;
} JOINT_6_DOF;

	//结构体类型的函数，返回一个结构体
	//solution取0就是关节角取正的解，反之取负的解
	//X-alpha,,Y-beta,,Z-gama
	JOINT_6_DOF inver_kinematic(float  x,float  y,float  z,float  alpha,float  beta,float  gama,uint8_t solution);	













#endif