#include "inver_kinematic.h"
#include "stm32f10x.h"
#include "IQmathLib.h"



#define r2d			57.2958
	

	
		
		

JOINT_6_DOF inver_kinematic(float  x,float  y,float  z,float  alpha,float  beta,float  gama,uint8_t solution)
{


		_iq iq_val_x,iq_val_y,iq_val_z,iq_val_cos,iq_val_r,iq_val_j1,iq_val_j2,iq_val_j3,iq_val_temp,iq_val_temp1,iq_val_temp2,iq_val_a,iq_val_b;
		
		_iq iq_val_r33,iq_val_r331,iq_val_r31,iq_val_r311,iq_val_r13,iq_val_r131,iq_val_r23,iq_val_r231,iq_val_r32,iq_val_r321,iq_val_z1,iq_val_y1,iq_val_z2,iq_val_z01,iq_val_y01,iq_val_z02;

		_iq iq_val_alpha,iq_val_beta,iq_val_gama,iq_val_n1,iq_val_n2,iq_val_n3,iq_val_h1,iq_val_h2,iq_val_h3,iq_val_o1,iq_val_o2,iq_val_o3;
		
		_iq last_y,last_z1,last_z2,Rzyz_13,Rzyz_23,Rzyz_33,wan_x,wan_y,wan_z;
	
	
	

	
	iq_val_alpha = _IQ(alpha);
	
	iq_val_beta = _IQ(beta);
	
	iq_val_gama = _IQ(gama); 
	
	JOINT_6_DOF joint_123;
	
	

	Rzyz_13 =_IQsin(  iq_val_beta );  

	Rzyz_23 =_IQmpy(-_IQsin( iq_val_alpha ),_IQcos( iq_val_beta )); 

	Rzyz_33 =_IQmpy( _IQcos( iq_val_beta ) ,_IQcos( iq_val_alpha )); 
	
	
	wan_x  = _IQmpy(Rzyz_13,_IQ(-4.2)); 
	wan_y  = _IQmpy(Rzyz_23,_IQ(-4.2)); 
	wan_z  = _IQmpy(Rzyz_33,_IQ(-4.2)); 
	

	iq_val_x  = wan_x + _IQ(x);	
	iq_val_y  = wan_y + _IQ(y);	
	iq_val_z  = wan_z + _IQ(z);	


	/*********************************/
	/**********求关节3的角度**********/
	/*********************************/

	//更改DH参数（改进的）时可以根据下面的公式计算
	//227.8 = a2^2+d3^2+L3^2
	//229.73 = 2*a2*L3
	iq_val_r = _IQmpy(iq_val_x,iq_val_x) + _IQmpy(iq_val_y,iq_val_y) + _IQmpy(iq_val_z,iq_val_z) - _IQ(229.73);  //tm1 =_IQtoF(iq_val_r);
	
	iq_val_cos = _IQdiv(iq_val_r,_IQ(227.8)); 
	
	if(solution == 0){iq_val_j3 = _IQacos(iq_val_cos);	}
	if(solution == 1){iq_val_j3 = -_IQacos(iq_val_cos);	}
	
	/*********************************/
	/**********求关节1的角度**********/
	/*********************************/

	iq_val_temp = _IQdiv(iq_val_y,iq_val_x);
	
	iq_val_j1 =  _IQatan(iq_val_temp );
	
	if(iq_val_y<0&iq_val_x==0){iq_val_j1 = -_IQatan(iq_val_temp );}
	
	/*********************************/
	/**********求关节2的角度**********/


	//计算公式：10是a2的长度,11.39是连杆3（L3）的长度
  
	iq_val_b = _IQ(-10.0) + _IQmpy(_IQ(-11.39) , _IQcos(iq_val_j3));
	
	iq_val_a = _IQmpy(_IQ(-11.39) , _IQsin(iq_val_j3));
	
	iq_val_temp = _IQmpy(iq_val_b,iq_val_b) + _IQmpy(iq_val_a,iq_val_a) -_IQmpy(iq_val_z,iq_val_z);
	
	iq_val_temp = iq_val_b +_IQsqrt(iq_val_temp);		
	
	iq_val_temp1 =  iq_val_a +iq_val_z;					
	
	if(iq_val_temp1==0){
		
		iq_val_temp = _IQdiv(iq_val_z,_IQ(38.8));
		
		iq_val_j2 = _IQatan(iq_val_temp);
	}
	else{
		iq_val_temp = _IQdiv(iq_val_temp,iq_val_temp1);
	
		iq_val_j2 = _IQatan(iq_val_temp);
	
		iq_val_j2 = _IQmpy(_IQ(2),iq_val_j2);
	}
	

	joint_123.jion_1 = _IQtoF(iq_val_j1)*r2d;
	
	joint_123.jion_2 = _IQtoF(iq_val_j2)*r2d;
	
	joint_123.jion_3 = _IQtoF(iq_val_j3)*r2d;
	
	
	
	
//下面是计算后三轴的转角
		
		//以下是期望姿态形成的矩阵移动角Rxyz
	iq_val_n1 = _IQmpy(_IQcos(iq_val_beta) , _IQcos(iq_val_gama));
	
	iq_val_n2 = _IQmpy(_IQcos(iq_val_gama),_IQmpy(_IQsin(iq_val_alpha),_IQsin(iq_val_beta)))+_IQmpy(_IQcos(iq_val_alpha),_IQsin(iq_val_gama));
	
	iq_val_n3 = -_IQmpy(_IQcos(iq_val_gama),_IQmpy(_IQcos(iq_val_alpha),_IQsin(iq_val_beta)))+_IQmpy(_IQsin(iq_val_alpha),_IQsin(iq_val_gama));
		
	iq_val_h1 = -_IQmpy(_IQcos(iq_val_beta),_IQsin(iq_val_gama));
	
	iq_val_h2 = -_IQmpy(_IQsin(iq_val_gama),_IQmpy(_IQsin(iq_val_alpha),_IQsin(iq_val_beta)))+_IQmpy(_IQcos(iq_val_alpha),_IQcos(iq_val_gama));
	
	iq_val_h3 = _IQmpy(_IQsin(iq_val_gama),_IQmpy(_IQcos(iq_val_alpha),_IQsin(iq_val_beta)))+_IQmpy(_IQsin(iq_val_alpha),_IQcos(iq_val_gama));
	
	iq_val_o1 = _IQsin(iq_val_beta);
	
	iq_val_o2 = -_IQmpy(_IQcos(iq_val_beta),_IQsin(iq_val_alpha));
	
	iq_val_o3 = _IQmpy(_IQcos(iq_val_beta),_IQcos(iq_val_alpha));



	//这是R03`的逆与期望姿态矩阵Rxyz的点乘
	iq_val_r31 = _IQmpy(iq_val_n1,_IQcos(iq_val_j1)) + _IQmpy(iq_val_n2,_IQsin(iq_val_j1));
	iq_val_r31 = _IQmpy(iq_val_r31,_IQcos(iq_val_j2+iq_val_j3)) ;
	iq_val_r31 = iq_val_r31 - _IQmpy(iq_val_n3,_IQsin(iq_val_j2+iq_val_j3)) ;
	
	iq_val_r33 = _IQmpy(iq_val_o1,_IQcos(iq_val_j1)) + _IQmpy(iq_val_o2,_IQsin(iq_val_j1));
	iq_val_r33 = _IQmpy(iq_val_r33,_IQcos(iq_val_j2+iq_val_j3)) ;
	iq_val_r33 = iq_val_r33 - _IQmpy(iq_val_o3,_IQsin(iq_val_j2+iq_val_j3)) ; 
	
	iq_val_r32 = _IQmpy(iq_val_h1,_IQcos(iq_val_j1)) + _IQmpy(iq_val_h2,_IQsin(iq_val_j1));
	iq_val_r32 = _IQmpy(iq_val_r32,_IQcos(iq_val_j2+iq_val_j3)) ;
	iq_val_r32 = iq_val_r32 - _IQmpy(iq_val_h3,_IQsin(iq_val_j2+iq_val_j3)) ;

	iq_val_r13 = _IQmpy(iq_val_o1,_IQcos(iq_val_j1)) + _IQmpy(iq_val_o2,_IQsin(iq_val_j1));
	iq_val_r13 = _IQmpy(iq_val_r13,_IQsin(iq_val_j2+iq_val_j3)) ;
	iq_val_r13 = -iq_val_r13 - _IQmpy(iq_val_o3,_IQcos(iq_val_j2+iq_val_j3)) ;
	
	iq_val_r23 = _IQmpy(iq_val_o2,_IQcos(iq_val_j1)) - _IQmpy(iq_val_o1,_IQsin(iq_val_j1));
	
//通过上述变换求解后三轴角度


/******************************末端的Y1轴********************************/
	if(_IQabs(iq_val_r33-_IQ(0))<_IQ(0.001)){iq_val_y1 =_IQ(1.57079);	}
	else {
		
	  iq_val_temp1 = _IQmag(iq_val_r31, iq_val_r32) ;	 
	  
		iq_val_temp1 = _IQdiv(iq_val_temp1,iq_val_r33); 
		
	  if(iq_val_r33<0){iq_val_y1 =_IQ(3.1415927)+_IQatan(iq_val_temp1);}
		
		
	  else{iq_val_y1 = _IQatan(iq_val_temp1);}
		

	}


	/*****************************末端的Z1轴************************************/
 if(_IQabs(iq_val_y1-_IQ(0))<_IQ(0.001)){ iq_val_z1=iq_val_gama;}		
 else {
		if(_IQabs(iq_val_y1-_IQ(3.1415927))<_IQ(0.001)){iq_val_z1=_IQ(0);}
		else{
			
	  iq_val_temp1 = _IQdiv(iq_val_r13,_IQsin(iq_val_y1)); 
	
 	  iq_val_temp2 = _IQdiv(iq_val_r23,_IQsin(iq_val_y1)); 
		
	  
			
		if(_IQabs(iq_val_temp1-_IQ(0))<_IQ(0.001)){
		if(iq_val_temp2<_IQ(0)){iq_val_z1 = -_IQ(1.5708);}
		else{iq_val_z1 = _IQ(1.5708);}
		}
		else{
						if(iq_val_temp1<0&&_IQabs(iq_val_temp2-_IQ(0))<_IQ(0.001))
				{									
//					iq_val_temp = _IQdiv(iq_val_temp2,iq_val_temp1);
//					iq_val_z1 = _IQ(3.1415926)-_IQatan(iq_val_temp);
					iq_val_z1 = _IQ(3.1415926);
				}
				else 
					{
					iq_val_temp = _IQdiv(iq_val_temp2,iq_val_temp1); 
					iq_val_z1 = _IQatan(iq_val_temp);
				 }
					
				 if(iq_val_temp1<_IQ(-0.003)&&iq_val_temp2>_IQ(0.001))
				 {
					 iq_val_temp = _IQdiv(iq_val_temp2,iq_val_temp1);
					iq_val_z1 = _IQ(3.1415926)+_IQatan(iq_val_temp);
					 
				 }
				 
				 if(iq_val_temp1<_IQ(-0.003)&&iq_val_temp2<_IQ(-0.002))
				 {
				  iq_val_temp = _IQdiv(iq_val_temp2,iq_val_temp1);
					iq_val_z1 = -_IQ(3.1415926)+_IQatan(iq_val_temp); 

				 }
			
		}
			
		 }}
	

		 
		 
		 /*****************************末端的Z2轴**************************************/
	if(_IQabs(iq_val_y1-_IQ(0))<_IQ(0.001)){ iq_val_z2=_IQ(0);}	
  else {
		if(_IQabs(iq_val_y1-_IQ(3.1415927))<_IQ(0.001)){iq_val_z2=iq_val_gama;}
		else{
			
	  iq_val_temp1 = _IQdiv(-iq_val_r31,_IQsin(iq_val_y1)); 
	
	  iq_val_temp2 = _IQdiv(iq_val_r32,_IQsin(iq_val_y1));
		
	  iq_val_temp = _IQdiv(iq_val_temp2,iq_val_temp1);  

			
				if(iq_val_temp1<=0&&_IQabs(iq_val_temp2-_IQ(0))<_IQ(0.001))
			{
//				iq_val_z2 = _IQ(3.1415926)-_IQatan(iq_val_temp);
				iq_val_z2 = _IQ(3.1415926);
			}else{
				iq_val_z2 = _IQatan(iq_val_temp);
			}
			
			if(iq_val_temp1<_IQ(-0.003)&&iq_val_temp2>_IQ(0.002))
				 {
					iq_val_z2 = _IQ(3.1415926)+_IQatan(iq_val_temp);
					 
				 }
				 
			if(iq_val_temp1<_IQ(-0.003)&&iq_val_temp2<_IQ(-0.002))
				 {
					iq_val_z2 = -_IQ(3.1415926)+_IQatan(iq_val_temp); 

				 }

		} }


	
	joint_123.ola_Y = _IQtoF(iq_val_y1)*r2d; 
	
	joint_123.ola_Z1 = _IQtoF(iq_val_z1)*r2d; 
	
	joint_123.ola_Z2 = _IQtoF(iq_val_z2)*r2d; 
	
	
	return ( joint_123 );
	
}








