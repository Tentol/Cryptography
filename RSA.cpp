#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath> 
#include<ctime>
#include<gmp.h>
#include<iostream>
using namespace std;

struct Key
{
	char* n;
	char* fn;
	char* p;
	char* q;
	char* d;
	int e;
};

Key *a;

void Square_and_Multiply(mpz_t z,const mpz_t x,const mpz_t c,const mpz_t n)
//z=x^c mod n
{
	char exp[2048];
	mpz_get_str(exp,2,c);
	mpz_t temp,pow;
	mpz_init(pow);
	mpz_init_set_ui(temp,1);
	mpz_mod(pow,x,n);
	for (int i=strlen(exp)-1;i>-1;i--)
	{
		if (exp[i]=='1')
		{
			mpz_mul(temp,temp,pow);
			mpz_mod(temp,temp,n);
		}
		mpz_mul(pow,pow,pow);
		mpz_mod(pow,pow,n);
	}
	mpz_set(z,temp);
	mpz_clear(temp);
	mpz_clear(pow);
}

int Miller_Rabin(mpz_t n, int k)
{
    mpz_t d, a, b, d_m;
    mpz_t pos1, pos2;
    mpz_init_set_ui(pos1,1);
    mpz_init_set(pos2,n);
    mpz_sub_ui(pos2,pos2,1);
    mpz_init(d);
    mpz_init(d_m);
    mpz_sub_ui(d, n, 1);
    mpz_sub_ui(d_m, n, 2);
    int r = 0;
    while (mpz_even_p(d)) 
	{
        mpz_cdiv_q_ui(d, d, 2);
        r++;
    }
    for (int i = 0; i < k; ++i) 
	{
        gmp_randstate_t grt;
        gmp_randinit_default(grt);
        gmp_randseed_ui(grt, time(NULL)+i);
        mpz_init(a);
        mpz_urandomb(a, grt, 1024);
        mpz_mod(a, a, d_m);
        mpz_add_ui(a, a, 2);
        mpz_init(b);
        mpz_powm(b, a, d, n); 
        int flag = 0;
        if (mpz_cmp(b,pos1) == 0 || mpz_cmp(b,pos2) == 0) 
			continue;
        for (int j = 0; j < r-1; ++j) 
		{
            mpz_mul(b,b,b);
            mpz_mod(b,b,n);
            if (mpz_cmp(b,pos1) == 0) return 0;
            if (mpz_cmp(b,pos2) == 0) 
			{
                flag=1;
                break;
            }
        }
        if (flag) continue;
        mpz_clear(d);
        mpz_clear(a);
        mpz_clear(b);
        mpz_clear(d_m);
        mpz_clear(pos1);
        mpz_clear(pos2);
        return 0;
    }
    mpz_clear(d);
    mpz_clear(a);
    mpz_clear(b);
	mpz_clear(d_m);
	mpz_clear(pos1);
	mpz_clear(pos2);
    return 1;
}

void invert(mpz_t &rop, mpz_t e, mpz_t f){
	mpz_t a, b, t1, t2, t, q;
	mpz_init_set(a, f);
	mpz_init_set(b, e);
	mpz_init(t1);             //t1=0 
	mpz_init_set_ui(t2, 1);   //t2=1
	mpz_init(t);
	mpz_init(q);
	
	while(mpz_cmp_ui(b, 0)!=0)
	{
		mpz_set(t, a);
		mpz_set(a, b);
		mpz_fdiv_qr(q,b,t,b);		// q=t/b  b=t%b
		
		mpz_mul(t,q,t2);              //t=q*t2
		mpz_sub(t,t1,t);        // t=t1-q*t2;
		mpz_set(t1, t2);
		mpz_set(t2, t); 
	}
	if(mpz_cmp_ui(t1,0)<0)
		mpz_add(t1, t1, f);
	mpz_set(rop, t1);
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(t1);
	mpz_clear(t2);
	mpz_clear(t);
	mpz_clear(q);
} 

mpz_t *Create_number()
{
	mpz_t n,m;
	int t=time(NULL);
	gmp_randstate_t state;
	mpz_init(n);
	mpz_init(m);
	gmp_randinit_default(state);
	gmp_randseed_ui(state,t);
	mpz_t *result=new mpz_t[2];
	mpz_init(result[0]);
	mpz_init(result[1]);
	mpz_urandomb(n,state,1024);
	if (mpz_even_p(n))
	{
		mpz_add_ui(n,n,1);
	}
	while (Miller_Rabin(n,10)!=1)
	{
		mpz_add_ui(n,n,1);
	}
	mpz_urandomb(m,state,1024);	
	if (mpz_even_p(m))
	{
		mpz_add_ui(m,m,1);
	}
	while (Miller_Rabin(m,10)!=1)
	{
		mpz_add_ui(m,m,1);
	}
	mpz_set(result[0],n);
	mpz_set(result[1],m);
	mpz_clear(n);
	mpz_clear(m);
	return result;
}

void ExtendEculid(mpz_t &Z, mpz_t a, mpz_t b, mpz_t &x, mpz_t &y, mpz_t c)
{
    mpz_t one,zero;
    mpz_init_set_ui(one,1);
    mpz_init_set_ui(zero,0);
    mpz_t temp1,temp2,result,mod;
    mpz_init(temp1);
    mpz_init(temp2);
    mpz_init(result);
    mpz_init(mod);

    mpz_t mul,div;
    mpz_init(mul);
    mpz_init(div);

    if(mpz_cmp(b,one) == 0)
     {
        mpz_set_ui(x,1);
        mpz_cdiv_q(temp1,c,a);
        mpz_ui_sub(y,0,temp1);
        mpz_set(Z, a);
        return;
     }
    mpz_mod(mod,a,b);
    ExtendEculid(result,b,mod,x,y,a);
    if(mpz_cmp(c,zero) != 0)
     {
        mpz_set(temp2,x);
        mpz_set(x,y);
        mpz_mul(mul,a,y);
        mpz_cdiv_q(div,c,mul);
        mpz_sub(y,temp2,div); 
     }
    mpz_set(Z, result);
    mpz_clear(one);
    mpz_clear(zero);
    return;
}

Key *Create_Key()
{
	mpz_t keyn,keye,keyf;
	mpz_t pp,qq;
	mpz_inits(keyn,keyf,pp,qq,0);
	mpz_init_set_ui(keye,65537);
	mpz_t *primes=Create_number();
	mpz_mul(keyn,primes[0],primes[1]); //n=pq;
	mpz_sub_ui(pp,primes[0],1);
	mpz_sub_ui(qq,primes[1],1);
	mpz_mul(keyf,pp,qq);//fn=(p-1)(q-1)
	mpz_t keyd;
	mpz_init(keyd);
	invert(keyd,keye,keyf);//求逆
	Key *result=new Key;
	char *buf_n=new char[1024];
	char *buf_f=new char[1024];
	char *buf_d=new char[1024];
	char *buf_p=new char[1024];
	char *buf_q=new char[1024];
	mpz_get_str(buf_n,10,keyn);
	result->n=buf_n;
	mpz_get_str(buf_f,10,keyf);
	result->fn=buf_f;
	mpz_get_str(buf_p,10,primes[0]);
	result->p=buf_p;
	mpz_get_str(buf_q,10,primes[1]);
	result->q=buf_q;
	mpz_get_str(buf_d,10,keyd);
	result->d=buf_d;
	result->e=65537;
	mpz_clear(keyn);
	mpz_clear(keye);
	mpz_clear(keyf);
	mpz_clear(pp);
	mpz_clear(qq);
	mpz_clear(keyd);
	return result;
}

char *encrypt_mod(const char *text,const char *keyn,int keye)
{
	mpz_t t,result,n,e;
	mpz_init_set_str(t,text,10);
	mpz_init_set_str(n,keyn,10);
	mpz_init_set_ui(result,0);
	mpz_init_set_ui(e,keye);
	Square_and_Multiply(result,t,e,n);
	char *r=new char[1025];
	mpz_get_str(r,10,result);
	mpz_clear(t);
	mpz_clear(result);
	mpz_clear(n);
	mpz_clear(e);
	return r;
}

char *decrypt_mod(const char *text,const char *keyn,const char *keyd)
{
	mpz_t t,result,n,d;
	mpz_init_set_str(t,text,10);
	mpz_init_set_str(n,keyn,10);
	mpz_init_set_str(d,keyd,10);
	mpz_init_set_ui(result,0);
	Square_and_Multiply(result,t,d,n);
	char *r=new char[1025];
	mpz_get_str(r,10,result);
	mpz_clear(t);
	mpz_clear(result);
	mpz_clear(n);
	mpz_clear(d);
	return r;
	
}

char *China_Remain(const char *P,const char *Q,const char *X,const char *D,const char *N)
{
	mpz_t M, C, p, q, d, e;
	mpz_t p_sub, q_sub, d_p, d_q, n;
	mpz_t m_p, m_q, p_f, q_f;
	mpz_init_set_str(C, X, 10);
	mpz_init_set_ui(e, 65537);
	mpz_init_set_str(n, N, 10);
	mpz_init_set_str(d, D, 10);
	mpz_init_set_str(p, P, 10);
	mpz_init_set_str(q, Q, 10);
	mpz_init(p_sub);
	mpz_init(q_sub);
	mpz_init(d_p);
	mpz_init(d_q);
	mpz_init(m_p);
	mpz_init(m_q);
	mpz_init(p_f);
	mpz_init(q_f);
	mpz_init(M);
	mpz_sub_ui(p_sub, p, 1);
	mpz_sub_ui(q_sub, q, 1);
	mpz_mod(d_p, d, p_sub);	
    mpz_mod(d_q, d, q_sub);
	mpz_powm(m_p, C, d_p, p);
	mpz_powm(m_q, C, d_q, q);
	invert(p_f, p, q);
	invert(q_f, q, p);
	mpz_mul(p, p, p_f);
	mpz_mul(q, q, q_f);
	mpz_mul(m_p, m_p, q);
	mpz_mul(m_q, m_q, p);
	mpz_add(M, m_p, m_q);
	mpz_mod(M, M, n);
	char *result=new char[1024];
	mpz_get_str(result,10, M);
	return result;
}

char *encrypt_China(const char *text,const char *keyp,const char *keyq,const char *keyn,int keye)
{
	mpz_t result;
	mpz_init_set_ui(result,0);
	char c[1025];
	_itoa(keye,c,10);
	char *r=new char[1025];
	r=China_Remain(keyp,keyq,text,c,keyn);
	return r;
}

char *decrypt_China(const char *text,const char *keyp,const char *keyq,const char *keyn,const char *keyd)
{
	mpz_t result;
	mpz_init_set_ui(result,0);
	char *r=new char[1024];
	r=China_Remain(keyp,keyq,text,keyd,keyn);
	return r;
} 

void MontMul(mpz_t result,const char *A,const char *B,const char *N)
{
	mpz_t a,b,n,tp,d;
	mpz_init_set_str(a,A,10);
	mpz_init_set_str(b,B,10);
	mpz_init_set_str(n,N,10);
	mpz_init_set_ui(tp,0);
	mpz_init_set_ui(d,1);
	while (mpz_cmp_ui(b,0))
	{
		if (mpz_odd_p(b))
		{
			mpz_mul(tp,d,a);
			mpz_mod(d,tp,n);
			mpz_sub_ui(b,b,1);
		}
		else
		{
			mpz_pow_ui(tp,a,2);
			mpz_mod(a,tp,n);
			mpz_divexact_ui(b,b,2);
		}
	}
	mpz_set(result,d);
}

char *encrypt_mont(const char *text,const char *keyn, int keye)
{
	mpz_t result;
	mpz_init(result);
	char c[2048];
	_itoa(keye,c,10);
	MontMul(result,text,c,keyn);
	char *r=new char[2048];
	mpz_get_str(r,10,result);
	return r;
}

char *decrypt_mont(const char *text,const char *keyn,const char *keyd)
{
	mpz_t result;
	mpz_init(result);
	MontMul(result,text,keyd,keyn);
	char *r=new char[2048];
	mpz_get_str(r,10,result);
	return r;
}

int main()
{
	int op=1;
	char *entxt,*detxt;
	while (op)
	{
		system("cls");
		printf("-------------------------\n");
		printf("1.密钥生成		2.模重复平方解密\n");
		printf("3.中国剩余定理解密	4.蒙哥马利解密\n");
		printf("0.退出\n");
		printf("请输入选择[0~4]:\n");
		cin>>op;
		switch (op)
		{
			case 1:
				{
					char text[2048];
					FILE *fp=fopen("1.txt","r");
					fgets(text,2048,fp);
					clock_t start,end;
					start=clock();
					a=Create_Key();
					end=clock();
					cout<<"密钥q=\n"<<a->p<<endl;
					cout<<"密钥p=\n"<<a->q<<endl;
					cout<<"n=\n"<<a->n<<endl;
					cout<<"φ(n)=\n"<<a->fn<<endl;
					cout<<"生成成功！用时："<<(double)end-start<<"ms"<<endl;
					entxt=encrypt_mod(text,a->n,a->e);
					cout<<"密文为："<<entxt<<endl; 
					fclose(fp);
					getchar();getchar();
					break;
				}
			case 2:
				{
					clock_t start=clock(); 
					for (int i=0;i<100;i++)
					{
						detxt=decrypt_mod(entxt,a->n,a->d);
					} 
					clock_t end=clock();
					cout<<"明文为："<<detxt<<endl;
					cout<<"解密完成！用时："<<((double)end-start)/100.00<<"ms"<<endl;
					getchar();getchar();
					break;
				}
			case 3:
				{
					clock_t start=clock(); 
					for (int i=0;i<100;i++)
					{
						char *detxt=decrypt_China(entxt,a->p,a->q,a->n,a->d);
					} 
					clock_t end=clock();
					cout<<"明文为："<<detxt<<endl;
					cout<<"解密完成！用时："<<((double)end-start)/100.00<<"ms"<<endl;
					getchar();getchar();
					break;
				}
			case 4:
				{
					clock_t start=clock();
					for (int i=0;i<100;i++)
					{
						char *detxt=decrypt_mont(entxt,a->n,a->d);
					} 
					clock_t end=clock();
					cout<<"明文为："<<detxt<<endl;
					cout<<"解密完成！用时："<<((double)end-start)/100.00<<"ms"<<endl;
					getchar();getchar();
					break;
				}
		
		} 
	} 
	//cout<<Create_Key()->d<<endl;
	
}
