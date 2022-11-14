#ifndef __INVER_KINEMATIC_H
#define __INVER_KINEMATIC_H 			   
#include "stm32f10x.h"

	
	typedef struct      //һ���ṹ��
{
  float jion_1;
  float jion_2;
	float jion_3;
	float ola_Z1;
  float ola_Y;
	float ola_Z2;
} JOINT_6_DOF;

	//�ṹ�����͵ĺ���������һ���ṹ��
	//solutionȡ0���ǹؽڽ�ȡ���Ľ⣬��֮ȡ���Ľ�
	//X-alpha,,Y-beta,,Z-gama
	JOINT_6_DOF inver_kinematic(float  x,float  y,float  z,float  alpha,float  beta,float  gama,uint8_t solution);	













#endif