#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#include <time.h>
#define pi atan(1)*4      
typedef struct{           
 double real;
 double imag;        //the structure to discribe plural
}plural;
            
void reverse(plural *x,int size_x);       //码位倒序               
void output(plural *x,int size_x);       //输出               
void add(plural,plural,plural *);    
void sub(plural,plural,plural *); 
void mul(plural,plural,plural *);   //复数运算
void fft(plural *x,int size_x);    //基二时域
void fft2(plural *x,int size_x);    //基二频域
void ifft(plural *x,int size_x);
void ifft2(plural *x,int size_x);
void linear_convolution(plural *,int,plural*,int,plural *);     //线性卷积
void fast_convolution1(plural *,int,plural*,int,plural *,int size_x);  // 基二时域的快速卷积
void fast_convolution2(plural *,int,plural*,int,plural *,int size_x);  // 基二频域的快速卷积
/*_^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^_*/


int main(void)               
{
 	int i;
   	plural x1[1024],y1[1024],c1[1024];
    plural x2[1024],y2[1024],c2[1024];
    plural x3[1024],y3[1024],c3[1024];
    int size_x=512;
    int  op1,end1,op2,end2,op3,end3;
    for(i=0;i<size_x;i++)
    { 
     x1[i].real=i+1;
     x2[i].real=i+1;
     x3[i].real=i+1;
     y1[i].real=i+3;
     y2[i].real=i+3;
     y3[i].real=i+3;
     x1[i].imag=i;
     y1[i].imag=i;
     x2[i].imag=i;
     y2[i].imag=i;
     x3[i].imag=i;
     y3[i].imag=i;
    }
            //以上是初始化序列


// 	printf("the answer of  fast_convolution2 is:\n");
//    op1=clock();
//    fast_convolution2(x1,2*size_x, y1,2*size_x, c1,2*size_x);
//    end1=clock();
//    output(c1,2*size_x);
    printf("\n\nthe answer of  fast_convolution1 is:\n");
    op2=clock();
    fast_convolution1(x2,2*size_x, y2,2*size_x, c2,2*size_x);
    end2=clock();
 output(c2,2*size_x);
 printf("\n\nthe answer for linear_convolution is : \n");
 op3=clock();
 linear_convolution(x3,size_x,y3,size_x,c3);   
 end3=clock();
 output(c3,2*size_x);
 printf("\n\n\n");
 printf("基2时=%dms，线性=%dms", end2-op2,end3-op3);
 } 


void fast_convolution1(plural *x,int m,plural *y,int n,plural *z,int size_x)
{
    int i;
    plural v;
 fft(x,m);
 fft(y,n);
 for(i=0;i<size_x;i++)     
 {
 mul(x[i],y[i],&v);
 z[i]=v;
 }
 ifft(z,m);
}
void fast_convolution2(plural *x,int m,plural *y,int n,plural *z,int size_x)
{
    int i;
    plural v;
 fft2(x,m);
 fft2(y,n);
 for(i=0;i<size_x;i++)
 {
 mul(x[i],y[i],&v);
 z[i]=v;
 }
 ifft2(z,m);
}






void reverse(plural *x,int size_x)                 //  码位倒序
{
 plural temp;
 unsigned short i=0,j=0,k=0;          
 double M;


 for(i=0;i<size_x;i++)
 {
 k=i;
 j=0;
 M=(log(size_x)/log(2));     // determine the series
 while((M--)>0)           //
 {
 j=j<<1;
 j|=(k & 1);      //如果k是奇数，那么j的最后一位 置1；
 k=k>>1;
 
 }
 if(j>i)
 {
 temp=x[i];
 x[i]=x[j];
 x[j]=temp;        //everything is ok;
 }
 }
}
 
 void output(plural *x,int size_x)            //序列输出                   
 {
   int i;
   for(i=0;i<size_x;i++)
   {
   printf("%.4f",x[i].real);  
        if(x[i].imag>=0)printf("+%.4fj  ",x[i].imag);  
        else //((x[i].imag)<0)printf("\n");  
        printf("%.4fj  ",x[i].imag);
  }
 printf("\n");
 }
 
 


void fft(plural *x,int size_x)     //基2时域 fft
{
 int i,j,m=0,n=0;
 plural A,B,C;
 plural *W; 
 reverse(x,size_x);
    W=(plural *)malloc(sizeof(plural) * size_x);
 for(i=0;i<size_x;i++)  
 {  
 W[i].real=cos(2*pi/size_x*i);   
 W[i].imag=-1*sin(2*pi/size_x*i);      //求  Wn() 旋转因子
 } 
 for(i=0;i<log(size_x)/log(2);i++)
 {
 n=1<<i;//or  pow(2,i)
 for(j=0;j<size_x;j=j+2*n)
 {
 for(m=0;m<n;m++)
 {
 mul(x[j+m+n],W[size_x*m/2/n],&C);
 add(x[j+m],C,&A);
 sub(x[j+m],C,&B);  
 x[j+m]=A;  
 x[j+m+n]=B; 
 }
 }
 }   
 
 
}


void fft2(plural *x,int size_x)      //基二频域fft
{
 int i,j,m=0,n=0;
 plural A,B,C;
 plural *W;
 W=(plural *)malloc(sizeof(plural) * size_x);
 for(i=0;i<size_x;i++)  
 {  
 W[i].real=cos(2*pi/size_x*i);   
 W[i].imag=-1*sin(2*pi/size_x*i);      //qiu  Wn()
 } 
 for(i=0;i<log(size_x)/log(2);i++)
 {
 n=size_x/2>>i;//or  pow(2,i)
 for(j=0;j<size_x;j=j+2*n)
 {
 for(m=0;m<n;m++)
 {
 add(x[j+m],x[j+m+n],&A);
 sub(x[j+m],x[j+m+n],&B);
 mul(B,W[size_x*m/2/n],&C); 
 x[j+m]=A;  
 x[j+m+n]=C; 
 }
 }
 } 
 reverse(x, size_x); 
 
}


void ifft(plural *x,int size_x)
{
 int i,j,m=0,n=0;
 plural A,B,C;
 plural *W; 
 plural k;
 k.real=1/size_x;
 k.imag=0;
 reverse(x,size_x);
    W=(plural *)malloc(sizeof(plural) * size_x);
 for(i=0;i<size_x;i++)  
 {  
 W[i].real=cos(2*pi/size_x*i);   
 W[i].imag=1*sin(2*pi/size_x*i);      //qiu  Wn()
 } 
 for(i=0;i<log(size_x)/log(2);i++)                //i 代表的是级
 {
 n=1<<i;//or  pow(2,i)  n是蝶形运算之间的偏移量
 for(j=0;j<size_x;j=j+2*n)                 
 {
 for(m=0;m<n;m++)
 {
 mul(x[j+m+n],W[size_x*m/2/n],&C);
 add(x[j+m],C,&A);
 sub(x[j+m],C,&B);
 x[j+m]=A;
 x[j+m+n]=B;
 }
 }
 } 
 for(i=0;i<size_x;i++)
 {
 x[i].imag=x[i].imag/size_x;
     x[i].real=x[i].real/size_x;
 }
        
}
void ifft2(plural *x,int size_x)
{
 int i,j,m=0,n=0;
 plural A,B,C;
 plural *W;
 W=(plural *)malloc(sizeof(plural) * size_x);
 for(i=0;i<size_x;i++)  
 {  
 W[i].real=cos(2*pi/size_x*i);   
 W[i].imag=1*sin(2*pi/size_x*i);      //qiu  Wn()
 } 
 for(i=0;i<log(size_x)/log(2);i++)
 {
 n=size_x/2>>i;//or  pow(2,i)
 for(j=0;j<size_x;j=j+2*n)
 {
 for(m=0;m<n;m++)
 {
 add(x[j+m],x[j+m+n],&A);
 sub(x[j+m],x[j+m+n],&B);
 mul(B,W[size_x*m/2/n],&C); 
 x[j+m]=A;  
 x[j+m+n]=C; 
 }
 }
 } 
 reverse(x, size_x);
 for(i=0;i<size_x;i++)
 {
 x[i].imag=x[i].imag/size_x;
     x[i].real=x[i].real/size_x;
 }
 
}






void linear_convolution(plural *x,int M,plural *h,int N,plural *z)    // 复数的线性卷积
{
 int n,m;
 plural A,B;
 B.real=0;
 B.imag=0;
 for(n=0;n<M+N;n++)
 {
 for(m=0;m<M;m++)
 {
 if(((n-m)>=0)&&((n-m)<N))
 {
 mul(x[m],h[n-m],&A);
 add(B, A, &B);
 }
 }
 z[n]=B;
 B.real=0;
 B.imag=0;
 }
 
}
    //复数的四则运算         


void add(plural a,plural b,plural *c)  
 {  
 c->real=a.real+b.real;
 c->imag=a.imag+b.imag;
 }
void sub(plural a,plural b,plural *c)
{
 c->real=a.real-b.real;
 c->imag=a.imag-b.imag;
 } 
void mul(plural a,plural b,plural *c)   
{
 c->real=a.real*b.real-a.imag*b.imag;
 c->imag=a.real*b.imag+a.imag*b.real;
 } 
