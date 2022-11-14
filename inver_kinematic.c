#include "inver_kinematic.h"
#include "stm32f10x.h"
#include "IQmathLib.h"



	 /*
	 *X_alpha：代表描述末端姿态的绕X轴的转角，是期望值，其他同理
	 *Joint_1：代表关节的转角，也就是舵机的转角
	 *cos_alpha：代表用期望转角计算的三角函数。其他同理
	 */

#define r2d			57.2958
	

		
		// float tm1,tm2,tm3,tm4,tm5,tm6;
		
		
	
//6轴还需要输入期望的姿态，这里设定未绕移动着XYZ的转角：alpha,beta gama
//但是此时不改了-8/18
JOINT_6_DOF inver_kinematic(float  x,float  y,float  z,float  alpha,float  beta,float  gama,uint8_t solution)
{

	//这些xyz的坐标是后三轴参考坐标系的位置坐标，并不是末端点的位置坐标,腕心位置
	//末端点位置需要经过计算才能得到腕心的位置，只有腕心的位置才能带入下面公式计算


		_iq iq_val_x,iq_val_y,iq_val_z,iq_val_cos,iq_val_r,iq_val_j1,iq_val_j2,iq_val_j3,iq_val_temp,iq_val_temp1,iq_val_temp2,iq_val_a,iq_val_b;
		
		_iq iq_val_r33,iq_val_r331,iq_val_r31,iq_val_r311,iq_val_r13,iq_val_r131,iq_val_r23,iq_val_r231,iq_val_r32,iq_val_r321,iq_val_z1,iq_val_y1,iq_val_z2,iq_val_z01,iq_val_y01,iq_val_z02;

		_iq iq_val_alpha,iq_val_beta,iq_val_gama,iq_val_n1,iq_val_n2,iq_val_n3,iq_val_h1,iq_val_h2,iq_val_h3,iq_val_o1,iq_val_o2,iq_val_o3;
		
		_iq last_y,last_z1,last_z2,Rzyz_13,Rzyz_23,Rzyz_33,wan_x,wan_y,wan_z;
	
	
	
	//y = y-5.5;
	
	iq_val_alpha = _IQ(alpha);
	
	iq_val_beta = _IQ(beta);
	
	iq_val_gama = _IQ(gama); 
	
	JOINT_6_DOF joint_123;
	
	

	
//下面的元素是Rzyz的，用来求腕心从原始位置经过旋转后的位置
/**10/11-算法优化处**/
	Rzyz_13 =_IQsin(  iq_val_beta );  

	Rzyz_23 =_IQmpy(-_IQsin( iq_val_alpha ),_IQcos( iq_val_beta )); //tm1 =_IQtoF(Rzyz_23);

	Rzyz_33 =_IQmpy( _IQcos( iq_val_beta ) ,_IQcos( iq_val_alpha )); 
	
	
	wan_x  = _IQmpy(Rzyz_13,_IQ(-4.2)); 
	wan_y  = _IQmpy(Rzyz_23,_IQ(-4.2)); 
	wan_z  = _IQmpy(Rzyz_33,_IQ(-4.2)); 
	
	//腕心旋转后的位置需要转换到机械臂的绝对坐标系中表示
	iq_val_x  = wan_x + _IQ(x);	//tm1 =_IQtoF(iq_val_x);
	iq_val_y  = wan_y + _IQ(y);	//tm1 =_IQtoF(iq_val_y);在第二代的机械臂中y坐标有一个偏移
	iq_val_z  = wan_z + _IQ(z);	//tm1 =_IQtoF(iq_val_z);
	//腕心的位置才是带入下面计算的参数

	/*********************************/
	/**********求关节3的角度**********/
	/*********************************/
	/*
	这角有两个解，从而导致逆运动学有两个解
	肘在上时（关节3取正）关节3的取值范围为（0，180），正正好时反cos的值域
	这个算法不存在分母为零的情况。一切都能正常运行
	*/
	//更改DH参数（改进的）时311和279.0也会发生变换，具体是多少需要在mathematical里面计算,也可以根据下面的公式计算
	//450 = a2^2+d3^2+L3^2
	//450 = 2*a2*L3
	iq_val_r = _IQmpy(iq_val_x,iq_val_x) + _IQmpy(iq_val_y,iq_val_y) + _IQmpy(iq_val_z,iq_val_z) - _IQ(229.73);  //tm1 =_IQtoF(iq_val_r);
	
	iq_val_cos = _IQdiv(iq_val_r,_IQ(227.8)); //tm2 =_IQtoF(iq_val_cos);
	
	if(solution == 0){iq_val_j3 = _IQacos(iq_val_cos);	}//取正号，就是肘关节在上的姿态
	if(solution == 1){iq_val_j3 = -_IQacos(iq_val_cos);	}//tm3 =_IQtoF(iq_val_j3);
	
	/*********************************/
	/**********求关节1的角度**********/
	/*********************************/
	/*
	从这个算法的角度看，关节1的取值范围是（-90，90）其实电机的转角和这个范围差不多
	只有一个特殊情况（奇异位置），当x=0时，分数就会变成无穷大按照经验看反正切后的值为90度
	但是这种情况下关节1就相当于固定角Rxyz种的z轴，这对于当前的姿态描述是无法处理的，所以应当
	避免这种情况。目前只能人为避免
	*/
	iq_val_temp = _IQdiv(iq_val_y,iq_val_x);
	
	iq_val_j1 =  _IQatan(iq_val_temp );
	
	if(iq_val_y<0&iq_val_x==0){iq_val_j1 = -_IQatan(iq_val_temp );}
	
	/*********************************/
	/**********求关节2的角度**********/
	/*********************************/
	/*
	逆运动学中本来就有些地方关节角度会极具变化
	这里的角度很迷，在关节三取正的情况下好像不会发生分母为零的情况但是
	会有误差极大的情况，因为a和b会变成很小的数但是他们的比值不小，
	这也许就是某些位置产生误差的原因吧？
	*/
	//-15和15也会随着DH参数而改变
	//计算公式：15是a2的长度,15是连杆3（L3）的长度
  
	iq_val_b = _IQ(-10.0) + _IQmpy(_IQ(-11.39) , _IQcos(iq_val_j3));
	
	iq_val_a = _IQmpy(_IQ(-11.39) , _IQsin(iq_val_j3));
	
	iq_val_temp = _IQmpy(iq_val_b,iq_val_b) + _IQmpy(iq_val_a,iq_val_a) -_IQmpy(iq_val_z,iq_val_z);
	
	iq_val_temp = iq_val_b +_IQsqrt(iq_val_temp);		//反正切函数参数的分子
	
	iq_val_temp1 =  iq_val_a +iq_val_z;					//反正切函数的分母
	
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
		/*这一部分是计算随前三轴末端位置而变化的姿态*/
		
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
	if(_IQabs(iq_val_r33-_IQ(0))<_IQ(0.001)){iq_val_y1 =_IQ(1.57079);	}//判断下式中的分母是否为零，特殊的r33为零式代表只绕z1轴旋转，这里只运动z1轴,且我规定y1轴只正转
	else {
		
	  iq_val_temp1 = _IQmag(iq_val_r31, iq_val_r32) ;	 //这里的值永远为正
	  
		iq_val_temp1 = _IQdiv(iq_val_temp1,iq_val_r33);  //tm2 =_IQtoF(iq_val_temp1);
		
	  if(iq_val_r33<0){iq_val_y1 =_IQ(3.1415927)+_IQatan(iq_val_temp1);}//这里需要四象限反正切函数因为y的取值范围为（0，180）
		//其实这并不是完整四象限反正切，只是x始终为正的情况，但正好满足y1轴的取值范围
		
	  else{iq_val_y1 = _IQatan(iq_val_temp1);}
		
//		tm1 =_IQtoF(iq_val_y1);
//		tm2 =_IQtoF(iq_val_r33);
//		tm3 =_IQtoF(iq_val_temp1);
	}




	/*****************************末端的Z1轴************************************/
	//sin(Y1)=两种情况，y1=0或者y1=180
 if(_IQabs(iq_val_y1-_IQ(0))<_IQ(0.001)){ iq_val_z1=iq_val_gama;}//zyz的y轴不转，只让z1轴转 		
 else {
		if(_IQabs(iq_val_y1-_IQ(3.1415927))<_IQ(0.001)){iq_val_z1=_IQ(0);}//这种情况很特殊z1不转动
		else{
			
	  iq_val_temp1 = _IQdiv(iq_val_r13,_IQsin(iq_val_y1)); //tm3 =_IQtoF(iq_val_temp1);
	
 	  iq_val_temp2 = _IQdiv(iq_val_r23,_IQsin(iq_val_y1)); //tm2 =_IQtoF(iq_val_temp2);
		
	  
			
		if(_IQabs(iq_val_temp1-_IQ(0))<_IQ(0.001)){
		if(iq_val_temp2<_IQ(0)){iq_val_z1 = -_IQ(1.5708);}//用来判断是90度还是-90度,只有在分母为零是进行判断，幸运的是分母为零时这个函数也可以正确运行
		else{iq_val_z1 = _IQ(1.5708);}
		}
		else{
						if(iq_val_temp1<0&&_IQabs(iq_val_temp2-_IQ(0))<_IQ(0.001))//虽然分母不为零但是当分子为零时整个分数也为零，符号就无法判断了
				{									//当分母小于零时处于第二象限
//					iq_val_temp = _IQdiv(iq_val_temp2,iq_val_temp1);
//					iq_val_z1 = _IQ(3.1415926)-_IQatan(iq_val_temp);
					iq_val_z1 = _IQ(3.1415926);
				}
				else 
					{
					iq_val_temp = _IQdiv(iq_val_temp2,iq_val_temp1); 
					iq_val_z1 = _IQatan(iq_val_temp);
				 }
					
				 if(iq_val_temp1<_IQ(-0.003)&&iq_val_temp2>_IQ(0.001))//在第二象限,这里不包括分子分母都为零的情况
				 {
					 iq_val_temp = _IQdiv(iq_val_temp2,iq_val_temp1);
					iq_val_z1 = _IQ(3.1415926)+_IQatan(iq_val_temp);
					 
				 }
				 
				 if(iq_val_temp1<_IQ(-0.003)&&iq_val_temp2<_IQ(-0.002))//在第三象限
				 {
				  iq_val_temp = _IQdiv(iq_val_temp2,iq_val_temp1);
					iq_val_z1 = -_IQ(3.1415926)+_IQatan(iq_val_temp); 

				 }
			
		}
			
//   
//		tm1 =_IQtoF(iq_val_temp1);
//		tm2 =_IQtoF(iq_val_temp2);
//		tm3 =_IQtoF(iq_val_y1);
//		tm4 =_IQtoF(iq_val_temp);
//		tm5 =_IQtoF(iq_val_z1);
		 }}
	

		 
		 
		 /*****************************末端的Z2轴**************************************/
	if(_IQabs(iq_val_y1-_IQ(0))<_IQ(0.001)){ iq_val_z2=_IQ(0);}	
  else {
		if(_IQabs(iq_val_y1-_IQ(3.1415927))<_IQ(0.001)){iq_val_z2=iq_val_gama;}//也许在这个值之间徘徊
		else{
			
	  iq_val_temp1 = _IQdiv(-iq_val_r31,_IQsin(iq_val_y1)); //tm1 =_IQtoF(iq_val_temp1);
	
	  iq_val_temp2 = _IQdiv(iq_val_r32,_IQsin(iq_val_y1)); //tm1 =_IQtoF(iq_val_temp2);
		
	  iq_val_temp = _IQdiv(iq_val_temp2,iq_val_temp1);   //tm1 =_IQtoF(iq_val_temp);

			
				if(iq_val_temp1<=0&&_IQabs(iq_val_temp2-_IQ(0))<_IQ(0.001))//位于第二象限,且分子为零
			{
//				iq_val_z2 = _IQ(3.1415926)-_IQatan(iq_val_temp);
				iq_val_z2 = _IQ(3.1415926);
			}else{
				iq_val_z2 = _IQatan(iq_val_temp);
			}
			
			if(iq_val_temp1<_IQ(-0.003)&&iq_val_temp2>_IQ(0.002))//在第二象限
				 {
					iq_val_z2 = _IQ(3.1415926)+_IQatan(iq_val_temp);
					 
				 }
				 
			if(iq_val_temp1<_IQ(-0.003)&&iq_val_temp2<_IQ(-0.002))//在第三象限
				 {
					iq_val_z2 = -_IQ(3.1415926)+_IQatan(iq_val_temp); 

				 }
			
			
	 // iq_val_z2 = _IQatan(iq_val_temp); //tm1 =_IQtoF(iq_val_z2);//这里应该不需要四象限反正切 
//		if(iq_val_temp1<_IQ(0)&&iq_val_temp2<_IQ(0)){iq_val_z2 = -(_IQ(3.141593)-iq_val_z2);}//第四象限
//		if(iq_val_temp1<_IQ(0)&&iq_val_temp2>_IQ(0)){iq_val_z2 = _IQ(3.141593)+iq_val_z2;}//第四象限
//		
//		if(iq_val_temp2<_IQ(0)&iq_val_temp1==_IQ(0)){iq_val_z2 = -iq_val_z2;}//用来判断是90度还是-90度,只有在分母为零是进行判断，幸运的是分母为零时这个函数也可以正确运行
//// 
//			
		//tm1 =_IQtoF(iq_val_z2);

		} }


	
	joint_123.ola_Y = _IQtoF(iq_val_y1)*r2d; 
	
	joint_123.ola_Z1 = _IQtoF(iq_val_z1)*r2d; 
	
	joint_123.ola_Z2 = _IQtoF(iq_val_z2)*r2d; 
	
	
	
	
	
	return ( joint_123 );
	
}








