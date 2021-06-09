//高通
#include <stdio.h>
#include <math.h>
#include <string.h>
#define fs 22000  //定义与课本相同的采样频率22000Hz 
#define pi 4.0*atan(1.0) 

int high_fliter(double h1[],double w[],double h[],double fH,double fL,int m); //h[n]就是滤波器 
int windows_hanning(double fh,double fl,double *w); //输入通带与阻带截止频率  输出w 
int windows_haming(double fh,double fl,double *w);
double div(double up,double down);


// 说明，在高通滤波器中呢，大部分与低通滤波器相同，稍微有修改，一下的是，high_fliter 与 low_fliter 略有所不同，
// 首先，利用高通的 通带截止频率与阻带截止频率 ，计算出相应阻带的频率， 在 在最后 多乘了一个 cos(n*omega0)罢了； 


int main(void)
{
 int i,N;
 double w[200]={0},h[200]={0},h1[200]={0},z[200]={0};
 printf("使用汉宁窗的滤波器序列h[n]为:\n");
 N=high_fliter(h1,w,h,8000,6000,0);   // 调用滤波函数，最后一个参数为 0 时，使用 汉宁窗 
 for(i=0;i<N;i++)
 {
 printf("%.4f ",h[i]);     // 打印滤波器序列 
  }
  printf("\n\n使用哈明窗的滤波器序列h[n]为:\n");
 memset(w,0,200*sizeof(double));     // 使用过的序列清零 
 memset(h1,0,200*sizeof(double)); 
 N=high_fliter(h1,w,z,8000,6000,1);   //
 for(i=0;i<N;i++)
 {
 printf("%.4f ",z[i]);
  } 
}
int windows_hanning(double fh,double fl,double *w) 
 {
  int TW,N;
  double x;
  int i,j;
  TW=fl-fh;                   
  x=3.32*fs/TW;      
 N=integer(x);         
 for(i=0;i<N;i++)
 {
 w[i]=0.5+0.5*cos(2*pi*(i-(N-1)/2)/(N-1));  
 }
 return N;            
 }
 
 int windows_haming(double fh,double fl,double *w)  // 哈明窗函数，与前面的汉宁窗基本一致，只是求的公式有些许不同而已。 
 {
  int TW,N;
  double x;
  int i,j;
  TW=fl-fh;
  x=3.44*fs/TW;
 N=integer(x);
 for(i=0;i<N;i++)
 {
 w[i]=0.54+0.46*cos(2*pi*(i-(N-1)/2)/(N-1));
 }
 return N;
 } 
 
 int high_fliter(double h1[],double w[],double h[],double fH,double fL,int m) // 滤波函数 
{
 int i=0,j,N,f,fH1,fL1;
 double f1,omega1;
 f=fs/2;
 fH1=f-fH;
 fL1=f-fL;
 f1=fH1+(fL1-fH1)/2;
 omega1=2*pi*f1/fs;     // 求出通带边缘截止数字频率 
 if(m==0)
 N=windows_hanning(fH1,fL1,w); // 满足m=0 时，调用 汉宁窗 
 else
 N=windows_haming(fH1,fL1,w);  // m!=0 时， 调用 哈明窗 
 printf("N=%d\n",N);
   for(i=0;i<N;i++)
   {
   // h1[i]=sin(omega1*(i-(N-1)/2))/((i-(N-1)/2)*pi);
     h1[i]=div(sin(omega1*(i-(N-1)/2)),(i-(N-1)/2)*pi); //求出 脉冲相应 
    h[i]=h1[i]*w[i]*cos(i*pi);            // 求出滤波器序列 
 }
   return N;
}


int integer(double x)      // 取序列长度为奇数的函数 
{
 int N;
 N=(int)x;
 if(N%2==0)
 {
 N=N+1;
 }
 else
 {
 N=N;  // 这里本应该是 N=N+2 才对， 但是考虑到应该与书上一致，故 写作 N=N； 
 }
 return N;
 }
 double div(double up,double down) // 排除了 sin0/0 的 除法函数 
 {
  double x;
  if(down==0)
  x=0.5;
  else 
  x=up/down;
  return x;
  
 } 



//低通

#include <stdio.h>
#include <math.h>
#include <string.h>
#define fs 10000  //定义与课本相同的采样频率10000Hz 
#define pi 4.0*atan(1.0) 
/********************************************************************************* 
******  low_fliter函数，h1作为滤波器脉冲响应，w 为窗函数 h 为输出的滤波器序列 
      fH fL 分别是通带截止频率，与阻带截止频率
  
  windows 函数， 也是输入两个截止频率 即可 fl为阻带 fh为通带  
***************************************************************************/ 
int low_fliter(double h1[],double w[],double h[],double fH,double fL,int m); //h[n]就是滤波器 
int windows_hanning(double fh,double fl,double *w); //输入通带与阻带截止频率  输出w 
int windows_haming(double fh,double fl,double *w);
double div(double up,double down); // 本函数是考虑到计算序列时， 会出现 sin0/0 的情况，在计算机下无法计算所想的解决办法。 

int main(void)
{
 int i,N;
 double x1[200]={0},x2[200]={0},x[200]={0};
 double f1=1000,f2=4500;    //定义两个函数的频率 
 double w[200]={0},h[200]={0},h1[200]={0},z[200]={0};
 for(i=0;i<20;i++)
 {
 x1[i]=sin(2*i*pi*f1/fs);
 x2[i]=sin(2*i*pi*f2/fs);
 x[i]=x1[i]+x2[i];        // 这里是产生两个不同频率的混合的序列，但在本实验中没有用到，在后期的进一步改造程序中会用到 
  }
  printf("使用汉宁窗的滤波器序列h[n]为:\n");
 N=low_fliter(h1,w,h,2000,3000,0);   // 调用滤波函数，最后一个参数为 0 时，使用 汉宁窗 
 for(i=0;i<N;i++)
 {
 printf("%.4f ",h[i]);     // 打印滤波器序列 
  } 
  printf("\n\n使用哈明窗的滤波器序列h[n]为:\n");
 memset(w,0,200*sizeof(double));     // 使用过的序列清零 
 memset(h1,0,200*sizeof(double)); 
 N=low_fliter(h1,w,z,2000,3000,1);   // 调用滤波函数，最后一个参数为 1 时，使用 哈明窗 
 for(i=0;i<N;i++)
 {
 printf("%.4f ",z[i]);
  } 
 
 } 
int windows_hanning(double fh,double fl,double *w) //汉宁窗函数 
 {
  int TW,N;
  double x;
  int i,j;
  TW=fl-fh;
  x=3.32*fs/TW;      //利用公式求出 序列长 x 为小数 
 N=integer(x);         // 调用函数取 序列长为奇数 
 for(i=0;i<N;i++)
 {
 w[i]=0.5+0.5*cos(2*pi*(i-(N-1)/2)/(N-1));   //根据相应的公式生成，窗函数 
 }
 return N;            //返回一个值，窗函数序列长度 
 }
int windows_rect(double fh,double fl,double *w) //汉宁窗函数 
 {
  int TW,N;
  double x;
  int i,j;
  TW=fl-fh;
  x=3.32*fs/TW;      //利用公式求出 序列长 x 为小数 
 N=integer(x);         // 调用函数取 序列长为奇数 
 for(i=0;i<N;i++)
 {
 w[i]=1;   //根据相应的公式生成，窗函数 
 }
 return N;            //返回一个值，窗函数序列长度 
 }
 
int windows_haming(double fh,double fl,double *w)  // 哈明窗函数，与前面的汉宁窗基本一致，只是求的公式有些许不同而已。 
 {
  int TW,N;
  double x;
  int i,j;
  TW=fl-fh;
  x=3.44*fs/TW;
 N=integer(x);
 for(i=0;i<N;i++)
 {
 w[i]=0.54+0.46*cos(2*pi*(i-(N-1)/2)/(N-1));
 }
 return N;
 } 
 


int low_fliter(double h1[],double w[],double h[],double fH,double fL,int m) // 滤波函数 
{
 int i=0,j,N;
 double f1,omega1;
 f1=fH+(fL-fH)/2;
 omega1=2*pi*f1/fs;     // 求出通带边缘截止数字频率 
 if(m==0)
 N=windows_hanning(fH,fL,w); // 满足m=0 时，调用 汉宁窗 
 else
 N=windows_haming(fH,fL,w);  // m!=0 时， 调用 哈明窗 
 printf("N=%d\n",N);
   for(i=0;i<N;i++)
   {
   // h1[i]=sin(omega1*(i-(N-1)/2))/((i-(N-1)/2)*pi);
     h1[i]=div(sin(omega1*(i-(N-1)/2)),(i-(N-1)/2)*pi); //求出 脉冲相应 
    h[i]=h1[i]*w[i];            // 求出滤波器序列 
 }
   return N;
}


 
int integer(double x)      // 取序列长度为奇数的函数 
{
 int N;
 N=(int)x;
 if(N%2==0)
 {
 N=N+1;
 }
 else
 {
 N=N;  // 这里本应该是 N=N+2 才对， 但是考虑到应该与书上一致，故 写作 N=N； 
 }
 return N;
 }
 double div(double up,double down) // 排除了 sin0/0 的 除法函数 
 {
  double x;
  if(down==0)
  x=0.5;
  else 
  x=up/down;
  return x;
  
 }