#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<iostream>
#include<fstream>
#include<time.h>
using namespace std; 

int Scase[16]={14,4,13,1,2,15,11,8,3,10,6,12,5,9,0,7};
int SScase[16]={14,4,13,8,2,15,11,3,10,0,6,12,5,9,7,1};
int deScase[16]={14,3,4,8,1,12,10,15,7,13,9,6,11,2,0,5};
int Pcase[16]={1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16};
int Key[5]={0x3A94,0xA94D,0x94D6,0x4D63,0xD63F};
int mkey[5];
int t1p[4],tp1[4];
int yy[128];
int SKey[128];
int SPcase[128]={1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64,67,70,73,76,79,82,85,88,91,94,97,100,103,106,109,112,115,118,121,124,
				2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92,95,98,101,104,107,110,113,116,119,122,125,
				3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102,105,108,111,114,117,120,123,126,128,127};
double time1,time2;
 
void Calculate_Key(int k)
{
	for (int i=0;i<=4;i++)
	{
		mkey[i]=((k<<(i*4))&0xFFFF0000)>>16;
	}
} 

int S(int x,int op)
{
	int a,b,c,d;
	a=(x&0xF000)>>12;
	b=(x&0x0F00)>>8;
	c=(x&0x00F0)>>4;
	d=x&0x000F;
	if (op==0)
	{
		a=Scase[a]<<12;
		b=Scase[b]<<8;
		c=Scase[c]<<4;
		d=Scase[d];
	}
	else if (op==1) 
	{
		a=deScase[a]<<12;
		b=deScase[b]<<8;
		c=deScase[c]<<4;
		d=deScase[d];
	}
	return a|b|c|d;
}

int _2to10(int *x)
{
	return (x[0]*8+x[1]*4+x[2]*2+x[3]);
}

void *_10to2(int b)
{
	t1p[3]=b%2;
	b=b>>1;
	t1p[2]=b%2;
	b=b>>1;
	t1p[1]=b%2;
	b=b>>1;
	t1p[0]=b%2;
}

void Super_SPN(int *x) 
{
	int w[128],u[128],v[128];
	int tp[4];
	int i,j;
	memset(SKey,0,sizeof(SKey));
	for (i=0;i<=127;i++)
		w[i]=x[i];
	for (int Nr=0;Nr<=8;Nr++)
	{
		for (i=0;i<=127;i++)
			u[i]=(w[i]+SKey[i])%2;
		for (i=0;i<=31;i++)
		{
			for (j=0;j<=3;j++)
				tp[j]=u[4*i+j];
			_10to2(Scase[_2to10(tp)]);
			for (j=0;j<=3;j++)
				v[4*i+j]=t1p[j];
		}
		for (i=0;i<=127;i++)
			w[i]=v[SPcase[i]];
	}
	for (i=0;i<=127;i++)
			u[i]=(w[i]+SKey[i])%2;
	for (i=0;i<=31;i++)
	{
		for (j=0;j<=3;j++)
			tp[j]=u[4*i+j];
		_10to2(Scase[_2to10(tp)]);
		for (j=0;j<=3;j++)
			v[4*i+j]=t1p[j];	
	}
	for (i=0;i<=127;i++)
	{
		w[i]=(v[i]+SKey[i])%2;
		yy[i]=w[i];
	}
}

int P(int x)
{
	int m[16],n[16];
	int c,y;
	c=x;
	y=0;
	for (int i=0;i<=15;i++)
	{
		m[15-i]=c%2;
		c=c>>1;
	}
	for (int i=0;i<=14;i++)
	{
		n[i]=m[Pcase[i]-1];
		y=y|n[i];
		y=y<<1;
	}
	n[15]=m[Pcase[15]-1];
	y=y|n[15];
	return y;
}

int BaseSpn(int x,int Nr,int *Key)
{
	int u,v,w,y;
	w=x;
	for (int i=1;i<=Nr-1;i++)
	{
		u=w^Key[i-1];
		v=S(u,0); 
		w=P(v);
 	}
 	u=w^Key[Nr-1];
 	v=S(u,0);
 	y=v^Key[Nr];
 	return y;
}

int deBaseSpn(int x,int Nr,int *Key)
{
	int u,v,w,y;
	v=x^Key[Nr];
	u=S(v,1);
	w=u^Key[Nr-1];
	for (int i=Nr-1;i>=1;i--)
	{
		v=P(w);
		u=S(v,1);
		w=u^Key[i-1];
	}
	return w;
}

int LinearAttack()
{
	//x为明文 y为密文 
	//K int Key[5]={0x3A94,0xA94D,0x94D6,0x4D63,0xD63F};
	int x,y,L1,L2;
	int y2,y4,x5,x7,x8,V2,V4,U2,U4,z;
	int u6,u8,u14,u16;
	int count[16][16];
	int mk[5]; 
	int flag=1,flag1,kk;
	int maxkey1,maxkey2,max=-1;
	int T=10000;
	memset(count,0,sizeof(count));
	kk=0;
	clock_t start=clock();
	for (int o=0;o<=T;o++)
	{
		x=rand()%65536;
		y=BaseSpn(x,4,Key);
		y2=(y&0x0F00)>>8;
		y4=y&0x000F;
		x5=(x&0x0800)>>11;
		x7=(x&0x0200)>>9;
		x8=(x&0x0100)>>8;
		for (L1=0;L1<=15;L1++)
		{
			for (L2=0;L2<=15;L2++)
			{
				V2=L1^y2;
				V4=L2^y4;
				U2=S(V2,1);
				U4=S(V4,1);
				u6=(U2&0x0004)>>2;
				u8=U2&0x0001;
				u14=(U4&0x0004)>>2;
				u16=U4&0x0001;
				z=x5^x7^x8^u6^u8^u14^u16;
				if (z==0)
				{
					count[L1][L2]++;
				}
			}
		}	
	}
	for (L1=0;L1<=15;L1++)
	{
		for (L2=0;L2<=15;L2++)
		{
			count[L1][L2]=abs(count[L1][L2]-T/2);
			if (count[L1][L2]>max)
			{
				max=count[L1][L2];
				maxkey1=L1;
				maxkey2=L2;
			}
		}
	}
	clock_t end=clock();
	time1=(double)(end-start);
	for (int k1=0;k1<=0xffff&&flag==1;k1++)
	{
		for (int k2=0;k2<=15&&flag==1;k2++)
		{
			for (int k3=0;k3<=15&&flag==1;k3++)
			{
				kk=(k1<<16)|(k2<<12)|(maxkey1<<8)|(k3<<4)|maxkey2;
				Calculate_Key(kk);
				int total=0; flag1=1;
				for (int turn=0;turn<10;turn++)
				{
					int plain=rand()%65535;
					int txt=BaseSpn(plain,4,Key);
					int text=BaseSpn(plain,4,mkey);
					if (txt!=text)
					{
						break;
					}
					if (turn==9)
						return kk;
				}
				
			}
		}
	}
}
void Create_SPcase()
{
	memset(SPcase,0,sizeof(SPcase));
	for (int i=0;i<=31;i++)
		SPcase[i]=1+4*i;
	for (int i=32;i<=63;i++)
		SPcase[i]=2+4*(i-32);
	for (int i=64;i<=95;i++)
		SPcase[i]=3+4*(i-64);
	for (int i=96;i<=127;i++)
		SPcase[i]=4+4*(i-96);
}

int PairAttack()
{
	int x,y,xx,yy,L1,L2;
	int xpair=0x0B00; 
	int mk[5]; 
	int flag=1,flag1,kk;
	int maxkey1,maxkey2,max=-1;
	int y1,yy1,y2,yy2,y3,yy3,y4,yy4;
	int v2,v4,u2,u4,vv2,vv4,uu2,uu4;
	int upair2,upair4;
	int count[16][16];
	int T=1000,num=0;
	memset(count,0,sizeof(count));
	clock_t start=clock();
	for (int k=0;k<=T;k++)
	{
		x=rand()%65536;
		y=BaseSpn(x,4,Key);
		xx=x^xpair;
		yy=BaseSpn(xx,4,Key);
		y1=(y&0xF000)>>12;
		yy1=(yy&0xF000)>>12;
		y2=(y&0x0F00)>>8;
		yy2=(yy&0x0F00)>>8;
		y3=(y&0x00F0)>>4;
		yy3=(yy&0x00F0)>>4;
		y4=y&0x000F;
		yy4=yy&0x000F; 
		if ((y1==yy1)&&(y3==yy3))
		{
			for (L1=0;L1<=15;L1++)
			{
				for (L2=0;L2<=15;L2++)
				{
					v2=L1^y2;
					v4=L2^y4;
					u2=S(v2,1);
					u4=S(v4,1);
					vv2=L1^yy2;
					vv4=L2^yy4;
					uu2=S(vv2,1);
					uu4=S(vv4,1);
					upair2=u2^uu2;
					upair4=u4^uu4;
					if ((upair2==6)&&(upair4==6))
					{
						count[L1][L2]++;
						num++;
					}
				}
			}
		}
	}
	for (L1=0;L1<=15;L1++)
	{
		for (L2=0;L2<=15;L2++)
		{
			if (count[L1][L2]>max)
			{
				max=count[L1][L2];
				maxkey1=L1;
				maxkey2=L2;
			}
		}
	}
	clock_t end=clock();
	time2=(double)(end-start);
	for (int k1=0;k1<=0xffff&&flag==1;k1++)
	{
		for (int k2=0;k2<=15&&flag==1;k2++)
		{
			for (int k3=0;k3<=15&&flag==1;k3++)
			{
				kk=(k1<<16)|(k2<<12)|(maxkey1<<8)|(k3<<4)|maxkey2;
				Calculate_Key(kk);
				int total=0; flag1=1;
				for (int turn=0;turn<10;turn++)
				{
					int plain=rand()%65535;
					int txt=BaseSpn(plain,4,Key);
					int text=BaseSpn(plain,4,mkey);
					if (txt!=text)
					{
						break;
					}
					if (turn==9)
						return kk;
				}
				
			}
		}
	}
}


int main()
{
	int op=1;
	while(op)
	{
		system("cls");
		printf("---------------------------\n");
		printf("1.SPN加密	2.SPN解密\n");
		printf("3.线性攻击	4.差分攻击\n"); 
		printf("5.SPN增强	0.退出\n");
		printf("---------------------------\n");
		printf("请输入选择[0~5]:\n");
		scanf("%d",&op);
		switch (op)
		{
			case 1:
				{
					printf("请以16进制输入要加密的明文：\n");
					int x;
					scanf("%x",&x); 
					printf("加密后的密文为：\n");
					printf("%x",BaseSpn(x,4,Key));
					getchar();getchar();
					break;
				}
			case 2:
				{
					printf("请以16进制输入密文：\n");
					int x;
					scanf("%x",&x);
					printf("加密前的明文为：\n");
					printf("%x",deBaseSpn(x,4,Key));
					getchar();getchar();
					break;
				}
			case 3:
				{
					clock_t start=clock();
					int result=LinearAttack();
					clock_t end=clock();
					printf("经线性分析，密钥可能为：%x\n",result); 
					printf("所用时间：%.lfms %.lfms\n",time1,(double)end-start);
					getchar();getchar();
					break;
				}
			case 4:
				{
					clock_t start=clock();
					int result=PairAttack();
					clock_t end=clock();
					printf("经差分分析，密钥可能为：%x\n",result); 
					printf("所用时间：%.lfms %.lfms\n",time2,(double)end-start);
					getchar();getchar();
					break;
				}
			case 5:
				{
					ifstream fin("2.txt");
					FILE *fp=fopen("D:\\3.txt","wb");
					int x[128],y[128];
					memset(yy,0,sizeof(yy));
					char r[16];
					char c;
					int k=0;
					clock_t start=clock();
					while(fin>>c)
					{
						x[k]=yy[k];
						k++;
						if (k==128)
						{
							Super_SPN(x);
							for (int i=0;i<16;i++)
							{
								r[i]=yy[i*8]*128+yy[i*8+1]*64+yy[i*8+2]*32+yy[i*8+3]*16+yy[i*8+4]*8+yy[i*8+5]*4+yy[i*8+6]*2+yy[i*8+7];
							}
							fwrite(r,1,16,fp);
							memset(x,0,sizeof(x));
							k=0;
						}
					}
					clock_t end=clock();
					fin.close();
					fclose(fp);
					printf("增强SPN加密成功!\n");
					getchar();getchar();
					break;
				}
			case 0:
				printf("欢迎下次再使用本系统！\n");
				break;
		}
	}
	return 0;
}
