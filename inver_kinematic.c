#include "inver_kinematic.h"
#include "stm32f10x.h"
#include "IQmathLib.h"



	 /*
	 *X_alpha����������ĩ����̬����X���ת�ǣ�������ֵ������ͬ��
	 *Joint_1������ؽڵ�ת�ǣ�Ҳ���Ƕ����ת��
	 *cos_alpha������������ת�Ǽ�������Ǻ���������ͬ��
	 */

#define r2d			57.2958
	

		
		// float tm1,tm2,tm3,tm4,tm5,tm6;
		
		
	
//6�ỹ��Ҫ������������̬�������趨δ���ƶ���XYZ��ת�ǣ�alpha,beta gama
//���Ǵ�ʱ������-8/18
JOINT_6_DOF inver_kinematic(float  x,float  y,float  z,float  alpha,float  beta,float  gama,uint8_t solution)
{

	//��Щxyz�������Ǻ�����ο�����ϵ��λ�����꣬������ĩ�˵��λ������,����λ��
	//ĩ�˵�λ����Ҫ����������ܵõ����ĵ�λ�ã�ֻ�����ĵ�λ�ò��ܴ������湫ʽ����


		_iq iq_val_x,iq_val_y,iq_val_z,iq_val_cos,iq_val_r,iq_val_j1,iq_val_j2,iq_val_j3,iq_val_temp,iq_val_temp1,iq_val_temp2,iq_val_a,iq_val_b;
		
		_iq iq_val_r33,iq_val_r331,iq_val_r31,iq_val_r311,iq_val_r13,iq_val_r131,iq_val_r23,iq_val_r231,iq_val_r32,iq_val_r321,iq_val_z1,iq_val_y1,iq_val_z2,iq_val_z01,iq_val_y01,iq_val_z02;

		_iq iq_val_alpha,iq_val_beta,iq_val_gama,iq_val_n1,iq_val_n2,iq_val_n3,iq_val_h1,iq_val_h2,iq_val_h3,iq_val_o1,iq_val_o2,iq_val_o3;
		
		_iq last_y,last_z1,last_z2,Rzyz_13,Rzyz_23,Rzyz_33,wan_x,wan_y,wan_z;
	
	
	
	//y = y-5.5;
	
	iq_val_alpha = _IQ(alpha);
	
	iq_val_beta = _IQ(beta);
	
	iq_val_gama = _IQ(gama); 
	
	JOINT_6_DOF joint_123;
	
	

	
//�����Ԫ����Rzyz�ģ����������Ĵ�ԭʼλ�þ�����ת���λ��
/**10/11-�㷨�Ż���**/
	Rzyz_13 =_IQsin(  iq_val_beta );  

	Rzyz_23 =_IQmpy(-_IQsin( iq_val_alpha ),_IQcos( iq_val_beta )); //tm1 =_IQtoF(Rzyz_23);

	Rzyz_33 =_IQmpy( _IQcos( iq_val_beta ) ,_IQcos( iq_val_alpha )); 
	
	
	wan_x  = _IQmpy(Rzyz_13,_IQ(-4.2)); 
	wan_y  = _IQmpy(Rzyz_23,_IQ(-4.2)); 
	wan_z  = _IQmpy(Rzyz_33,_IQ(-4.2)); 
	
	//������ת���λ����Ҫת������е�۵ľ�������ϵ�б�ʾ
	iq_val_x  = wan_x + _IQ(x);	//tm1 =_IQtoF(iq_val_x);
	iq_val_y  = wan_y + _IQ(y);	//tm1 =_IQtoF(iq_val_y);�ڵڶ����Ļ�е����y������һ��ƫ��
	iq_val_z  = wan_z + _IQ(z);	//tm1 =_IQtoF(iq_val_z);
	//���ĵ�λ�ò��Ǵ����������Ĳ���

	/*********************************/
	/**********��ؽ�3�ĽǶ�**********/
	/*********************************/
	/*
	����������⣬�Ӷ��������˶�ѧ��������
	������ʱ���ؽ�3ȡ�����ؽ�3��ȡֵ��ΧΪ��0��180����������ʱ��cos��ֵ��
	����㷨�����ڷ�ĸΪ��������һ�ж�����������
	*/
	//����DH�������Ľ��ģ�ʱ311��279.0Ҳ�ᷢ���任�������Ƕ�����Ҫ��mathematical�������,Ҳ���Ը�������Ĺ�ʽ����
	//450 = a2^2+d3^2+L3^2
	//450 = 2*a2*L3
	iq_val_r = _IQmpy(iq_val_x,iq_val_x) + _IQmpy(iq_val_y,iq_val_y) + _IQmpy(iq_val_z,iq_val_z) - _IQ(229.73);  //tm1 =_IQtoF(iq_val_r);
	
	iq_val_cos = _IQdiv(iq_val_r,_IQ(227.8)); //tm2 =_IQtoF(iq_val_cos);
	
	if(solution == 0){iq_val_j3 = _IQacos(iq_val_cos);	}//ȡ���ţ�������ؽ����ϵ���̬
	if(solution == 1){iq_val_j3 = -_IQacos(iq_val_cos);	}//tm3 =_IQtoF(iq_val_j3);
	
	/*********************************/
	/**********��ؽ�1�ĽǶ�**********/
	/*********************************/
	/*
	������㷨�ĽǶȿ����ؽ�1��ȡֵ��Χ�ǣ�-90��90����ʵ�����ת�Ǻ������Χ���
	ֻ��һ���������������λ�ã�����x=0ʱ�������ͻ���������վ��鿴�����к��ֵΪ90��
	������������¹ؽ�1���൱�ڹ̶���Rxyz�ֵ�z�ᣬ����ڵ�ǰ����̬�������޷�����ģ�����Ӧ��
	�������������Ŀǰֻ����Ϊ����
	*/
	iq_val_temp = _IQdiv(iq_val_y,iq_val_x);
	
	iq_val_j1 =  _IQatan(iq_val_temp );
	
	if(iq_val_y<0&iq_val_x==0){iq_val_j1 = -_IQatan(iq_val_temp );}
	
	/*********************************/
	/**********��ؽ�2�ĽǶ�**********/
	/*********************************/
	/*
	���˶�ѧ�б�������Щ�ط��ؽڽǶȻἫ�߱仯
	����ĽǶȺ��ԣ��ڹؽ���ȡ��������º��񲻻ᷢ����ĸΪ����������
	����������������Ϊa��b���ɺ�С�����������ǵı�ֵ��С��
	��Ҳ�����ĳЩλ�ò�������ԭ��ɣ�
	*/
	//-15��15Ҳ������DH�������ı�
	//���㹫ʽ��15��a2�ĳ���,15������3��L3���ĳ���
  
	iq_val_b = _IQ(-10.0) + _IQmpy(_IQ(-11.39) , _IQcos(iq_val_j3));
	
	iq_val_a = _IQmpy(_IQ(-11.39) , _IQsin(iq_val_j3));
	
	iq_val_temp = _IQmpy(iq_val_b,iq_val_b) + _IQmpy(iq_val_a,iq_val_a) -_IQmpy(iq_val_z,iq_val_z);
	
	iq_val_temp = iq_val_b +_IQsqrt(iq_val_temp);		//�����к��������ķ���
	
	iq_val_temp1 =  iq_val_a +iq_val_z;					//�����к����ķ�ĸ
	
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
	
	
	
	
	
	
	
	
	
	
	
	
//�����Ǽ���������ת��
		/*��һ�����Ǽ�����ǰ����ĩ��λ�ö��仯����̬*/
		
		//������������̬�γɵľ����ƶ���Rxyz
	iq_val_n1 = _IQmpy(_IQcos(iq_val_beta) , _IQcos(iq_val_gama));
	
	iq_val_n2 = _IQmpy(_IQcos(iq_val_gama),_IQmpy(_IQsin(iq_val_alpha),_IQsin(iq_val_beta)))+_IQmpy(_IQcos(iq_val_alpha),_IQsin(iq_val_gama));
	
	iq_val_n3 = -_IQmpy(_IQcos(iq_val_gama),_IQmpy(_IQcos(iq_val_alpha),_IQsin(iq_val_beta)))+_IQmpy(_IQsin(iq_val_alpha),_IQsin(iq_val_gama));
		
	iq_val_h1 = -_IQmpy(_IQcos(iq_val_beta),_IQsin(iq_val_gama));
	
	iq_val_h2 = -_IQmpy(_IQsin(iq_val_gama),_IQmpy(_IQsin(iq_val_alpha),_IQsin(iq_val_beta)))+_IQmpy(_IQcos(iq_val_alpha),_IQcos(iq_val_gama));
	
	iq_val_h3 = _IQmpy(_IQsin(iq_val_gama),_IQmpy(_IQcos(iq_val_alpha),_IQsin(iq_val_beta)))+_IQmpy(_IQsin(iq_val_alpha),_IQcos(iq_val_gama));
	
	iq_val_o1 = _IQsin(iq_val_beta);
	
	iq_val_o2 = -_IQmpy(_IQcos(iq_val_beta),_IQsin(iq_val_alpha));
	
	iq_val_o3 = _IQmpy(_IQcos(iq_val_beta),_IQcos(iq_val_alpha));



	//����R03`������������̬����Rxyz�ĵ��
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
	
//ͨ�������任��������Ƕ�


/******************************ĩ�˵�Y1��********************************/
	if(_IQabs(iq_val_r33-_IQ(0))<_IQ(0.001)){iq_val_y1 =_IQ(1.57079);	}//�ж���ʽ�еķ�ĸ�Ƿ�Ϊ�㣬�����r33Ϊ��ʽ����ֻ��z1����ת������ֻ�˶�z1��,���ҹ涨y1��ֻ��ת
	else {
		
	  iq_val_temp1 = _IQmag(iq_val_r31, iq_val_r32) ;	 //�����ֵ��ԶΪ��
	  
		iq_val_temp1 = _IQdiv(iq_val_temp1,iq_val_r33);  //tm2 =_IQtoF(iq_val_temp1);
		
	  if(iq_val_r33<0){iq_val_y1 =_IQ(3.1415927)+_IQatan(iq_val_temp1);}//������Ҫ�����޷����к�����Ϊy��ȡֵ��ΧΪ��0��180��
		//��ʵ�Ⲣ�������������޷����У�ֻ��xʼ��Ϊ�������������������y1���ȡֵ��Χ
		
	  else{iq_val_y1 = _IQatan(iq_val_temp1);}
		
//		tm1 =_IQtoF(iq_val_y1);
//		tm2 =_IQtoF(iq_val_r33);
//		tm3 =_IQtoF(iq_val_temp1);
	}




	/*****************************ĩ�˵�Z1��************************************/
	//sin(Y1)=���������y1=0����y1=180
 if(_IQabs(iq_val_y1-_IQ(0))<_IQ(0.001)){ iq_val_z1=iq_val_gama;}//zyz��y�᲻ת��ֻ��z1��ת 		
 else {
		if(_IQabs(iq_val_y1-_IQ(3.1415927))<_IQ(0.001)){iq_val_z1=_IQ(0);}//�������������z1��ת��
		else{
			
	  iq_val_temp1 = _IQdiv(iq_val_r13,_IQsin(iq_val_y1)); //tm3 =_IQtoF(iq_val_temp1);
	
 	  iq_val_temp2 = _IQdiv(iq_val_r23,_IQsin(iq_val_y1)); //tm2 =_IQtoF(iq_val_temp2);
		
	  
			
		if(_IQabs(iq_val_temp1-_IQ(0))<_IQ(0.001)){
		if(iq_val_temp2<_IQ(0)){iq_val_z1 = -_IQ(1.5708);}//�����ж���90�Ȼ���-90��,ֻ���ڷ�ĸΪ���ǽ����жϣ����˵��Ƿ�ĸΪ��ʱ�������Ҳ������ȷ����
		else{iq_val_z1 = _IQ(1.5708);}
		}
		else{
						if(iq_val_temp1<0&&_IQabs(iq_val_temp2-_IQ(0))<_IQ(0.001))//��Ȼ��ĸ��Ϊ�㵫�ǵ�����Ϊ��ʱ��������ҲΪ�㣬���ž��޷��ж���
				{									//����ĸС����ʱ���ڵڶ�����
//					iq_val_temp = _IQdiv(iq_val_temp2,iq_val_temp1);
//					iq_val_z1 = _IQ(3.1415926)-_IQatan(iq_val_temp);
					iq_val_z1 = _IQ(3.1415926);
				}
				else 
					{
					iq_val_temp = _IQdiv(iq_val_temp2,iq_val_temp1); 
					iq_val_z1 = _IQatan(iq_val_temp);
				 }
					
				 if(iq_val_temp1<_IQ(-0.003)&&iq_val_temp2>_IQ(0.001))//�ڵڶ�����,���ﲻ�������ӷ�ĸ��Ϊ������
				 {
					 iq_val_temp = _IQdiv(iq_val_temp2,iq_val_temp1);
					iq_val_z1 = _IQ(3.1415926)+_IQatan(iq_val_temp);
					 
				 }
				 
				 if(iq_val_temp1<_IQ(-0.003)&&iq_val_temp2<_IQ(-0.002))//�ڵ�������
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
	

		 
		 
		 /*****************************ĩ�˵�Z2��**************************************/
	if(_IQabs(iq_val_y1-_IQ(0))<_IQ(0.001)){ iq_val_z2=_IQ(0);}	
  else {
		if(_IQabs(iq_val_y1-_IQ(3.1415927))<_IQ(0.001)){iq_val_z2=iq_val_gama;}//Ҳ�������ֵ֮���ǻ�
		else{
			
	  iq_val_temp1 = _IQdiv(-iq_val_r31,_IQsin(iq_val_y1)); //tm1 =_IQtoF(iq_val_temp1);
	
	  iq_val_temp2 = _IQdiv(iq_val_r32,_IQsin(iq_val_y1)); //tm1 =_IQtoF(iq_val_temp2);
		
	  iq_val_temp = _IQdiv(iq_val_temp2,iq_val_temp1);   //tm1 =_IQtoF(iq_val_temp);

			
				if(iq_val_temp1<=0&&_IQabs(iq_val_temp2-_IQ(0))<_IQ(0.001))//λ�ڵڶ�����,�ҷ���Ϊ��
			{
//				iq_val_z2 = _IQ(3.1415926)-_IQatan(iq_val_temp);
				iq_val_z2 = _IQ(3.1415926);
			}else{
				iq_val_z2 = _IQatan(iq_val_temp);
			}
			
			if(iq_val_temp1<_IQ(-0.003)&&iq_val_temp2>_IQ(0.002))//�ڵڶ�����
				 {
					iq_val_z2 = _IQ(3.1415926)+_IQatan(iq_val_temp);
					 
				 }
				 
			if(iq_val_temp1<_IQ(-0.003)&&iq_val_temp2<_IQ(-0.002))//�ڵ�������
				 {
					iq_val_z2 = -_IQ(3.1415926)+_IQatan(iq_val_temp); 

				 }
			
			
	 // iq_val_z2 = _IQatan(iq_val_temp); //tm1 =_IQtoF(iq_val_z2);//����Ӧ�ò���Ҫ�����޷����� 
//		if(iq_val_temp1<_IQ(0)&&iq_val_temp2<_IQ(0)){iq_val_z2 = -(_IQ(3.141593)-iq_val_z2);}//��������
//		if(iq_val_temp1<_IQ(0)&&iq_val_temp2>_IQ(0)){iq_val_z2 = _IQ(3.141593)+iq_val_z2;}//��������
//		
//		if(iq_val_temp2<_IQ(0)&iq_val_temp1==_IQ(0)){iq_val_z2 = -iq_val_z2;}//�����ж���90�Ȼ���-90��,ֻ���ڷ�ĸΪ���ǽ����жϣ����˵��Ƿ�ĸΪ��ʱ�������Ҳ������ȷ����
//// 
//			
		//tm1 =_IQtoF(iq_val_z2);

		} }


	
	joint_123.ola_Y = _IQtoF(iq_val_y1)*r2d; 
	
	joint_123.ola_Z1 = _IQtoF(iq_val_z1)*r2d; 
	
	joint_123.ola_Z2 = _IQtoF(iq_val_z2)*r2d; 
	
	
	
	
	
	return ( joint_123 );
	
}








