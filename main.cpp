#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <sstream>

using namespace std;

const int W=800;
const int H=600;
const double e=2.7182818284; //liczba Eulera
const double pi=3.1415926535;

int N=300;
int size_=3000;

float L1 = W/12, L2=W/12, m1 = 1, m2 = 1;
float px0 = W/2, py0 = H/4;
float dt = 0.1;
float g=100;

int R[W][H];
int G[W][H];
int B[W][H];

void saveppm(string name, string kom=""); // eksport do pliku ppm
void clc();
double exponent(float k, int phi,float half);
void Losowe();
void Rysunek_Prep(float a, float k,int scale);
void Rysunek(float a, float k,int scale);
void line(int x0, int y0, int xk, int yk, float grain, int r, int g, int b);
void Animacja();

int main()
{


//    Losowe();
//    saveppm("Lista9_Zadanie1.ppm", "Kolory pikseli generowane pseudolosowo");

//    Rysunek_Prep(1,0.3,5);
//    Rysunek(1,0.3,5);
//    saveppm("Lista9_Zadanie2.ppm");

    Animacja();

    return 0;
}

void saveppm(string name, string kom) // eksport do pliku ppm
{
    ofstream file(name);

    file << "P3" << endl << "#" << kom << endl << W << " " << H << endl << 255 << endl;

    for(int j=0;j<H;j++)
        {
            for(int i=0;i<W;i++)
                {
                    file << R[i][j] << " " << G[i][j] << " " << B[i][j] << " ";
                }
                file << endl;
        }
    file.close();
}
void clc()
{
    for(int i=0;i<W;i++)
        for(int j=0;j<H;j++)
             {
                 R[i][j] = G[i][j] = B[i][j] = 0;
             }
}
double exponent(float k, int phi,float half)
{
    return 1 + k*phi/half + (k*phi/half)*(k*phi/half)/2 + (k*phi/half)*(k*phi/half)*(k*phi/half)/6 + (k*phi/half)*(k*phi/half)*(k*phi/half)*(k*phi/half)/24;// + (k*phi/half)*(k*phi/half)*(k*phi/half)*(k*phi/half)*(k*phi/half)/120;
}
void Losowe()
{
    clc();
    for(int i=0;i<W;i++)
        for(int j=0;j<H;j++)
             {
                 srand(i*j*time(NULL));
                 R[i][j] = int(float(rand())/RAND_MAX*255.0);
                 G[i][j] = int(float(rand())/RAND_MAX*255.0);
                 B[i][j] = int(float(rand())/RAND_MAX*255.0);
             }
}
void Rysunek_Prep(float a, float k,int scale)
{
    clc();
    Losowe();

    static int steps = 50000;
    static float half = steps/2.;

    static float small_step=0.5/255.;

    int x,y;


    for(int phi=0; phi<=steps*30;phi++)
    {

        for(int i=0;i<100;i++)
            {
            x=int(-2*scale*(a-i*small_step)*cos(-phi/half*pi)*exponent(-k,phi*pi,half) + W/2);
            y=int(-2*scale*(a-i*small_step)*sin(-phi/half*pi)*exponent(-k,phi*pi,half) + H/2);

            if(x>=0&&x<=W&&y>=0&&y<=H)
                {
                    R[x][y] = 0;
                    G[x][y] = 0;
                    B[x][y] = 0;
                }
            }
    }
}
void Rysunek(float a, float k,int scale)
{
    static int steps = 50000;
    static float half = steps/2.;

    static float small_step=0.5/255.;

    int x,y;

    cout<< "generuje "<<a<<":"<<k<<endl;



    for(int phi=0; phi<=steps*30;phi++)
    {

        for(int i=0;i<100;i++)
            {
            x=int(scale*(a+i*small_step)*cos(phi/half*pi)*exponent(k,phi*pi,half)) + int(W/2);
            y=int(scale*(a+i*small_step)*sin(phi/half*pi)*exponent(k,phi*pi,half)) + int(H/2);

            if(x>=0&&x<=W&&y>=0&&y<=H)
                {
                    R[x][y] = 255-int(255*float(i/100.));
                    G[x][y] = int(float(i/100.)*255*2*abs(x/W-W/2));//-int(i/2.);
                    B[x][y] = int(float(i/100.)*255*2*abs(y/H-H/2));
                }
            }
    }
    cout<< "narysowano "<<a<<":"<<k<<endl<<endl;

}
void line(int x0, int y0, int xk, int yk,float grain, int r, int g, int b)
{
    double a=float(yk-y0)/(xk-x0);
    double b0=double(y0)-a*x0;

    int x, y;
    if(abs(xk-x0)>=0.05)
        {
            for(int i=0;i<=abs(xk-x0)*grain;i++)
            {
                int sign;
                if(xk>=x0)
                    sign=1;
                else
                    sign=-1;
                x = int(x0+sign*i/grain);
                y = int(a*(x0+sign*i/grain)+b0);

                if(x>=0&&x<=W&&y>=0&&y<=H)
                    {
                        R[x][y] = r;
                        G[x][y] = g;
                        B[x][y] = b;
                    }
            }
        }
        else
            for(int i=0;i<=abs(yk-y0)*grain;i++)
            {

                x = x0;
                y = int(y0+i/grain);

                if(x>=0&&x<=W&&y>=0&&y<=H)
                    {
                        R[x][y] = r;
                        G[x][y] = g;
                        B[x][y] = b;
                    }
            }
}
void Animacja()
{
    clc();

    int x1, y1;
    int x2, y2;

    double theta1[N];
    double theta2[N];

    int r,gr,b;

    double dtheta1[N], ddtheta1[N];
    double dtheta2[N], ddtheta2[N];


    for(int n=0;n<N;n++)
    {
        theta1[n] = pi/4 + 0.00001*double(n);
        theta2[n] = pi/3 + 0.00001*double(n);

        cout<<theta1[n]<<"\n";

        dtheta1[n] = 0.00;
        ddtheta1[n] = 0.00;

        dtheta2[n] = 0.00;
        ddtheta2[n] = 0.00;
    }

    for(int i=0; i<=size_; i++)
    {
        clc();

        for(int n=0;n<N;n++)
        {
            //cout<<n<<"\'"<<theta1[n]<<":"<<dtheta1[n]<<".."<<ddtheta1[n]<<"\n";
            //calculating angles

            theta1[n]=theta1[n] + dtheta1[n]*dt;
            while(theta1[n]>2*pi)
                {theta1[n]=theta1[n]-2*pi;}
            while(theta1[n]<-2*pi)
                {theta1[n]=theta1[n]+2*pi;}

            theta2[n]=theta2[n] + dtheta2[n]*dt;
            while(theta2[n]>2*pi)
                {theta2[n]=theta2[n]-2*pi;}
            while(theta2[n]<-2*pi)
                {theta2[n]=theta2[n]+2*pi;}

            ddtheta1[n]=(-g*(2*m1 + m2)*sin(theta1[n]) - m2*g*sin(theta1[n] - 2*theta2[n])
                         - 2*sin(theta1[n] - theta2[n])*m2*(dtheta2[n]*dtheta2[n]*L2 + dtheta1[n]*dtheta1[n]*L1*cos(theta1[n] - theta2[n])));
            ddtheta1[n]/=L1*(2*m1 + m2 - m2*cos(2*theta1[n] - 2*theta2[n]));

            while(ddtheta1[n]>=2*pi/(dt*dt))
                ddtheta1[n]=ddtheta1[n]-2*pi/(dt*dt);
            while(ddtheta1[n]<=-2*pi/(dt*dt))
                ddtheta1[n]=ddtheta1[n]+2*pi/(dt*dt);

            ddtheta2[n]=(2*sin(theta1[n] - theta2[n])*(dtheta1[n]*dtheta1[n]*L1*(m1 + m2)+g*(m1 + m2)*cos(theta1[n]) + dtheta2[n]*dtheta2[n]*L2*m2*cos(theta1[n] - theta2[n])));
            ddtheta2[n]/=L2*(2*m1 + m2 - m2*cos(2*theta1[n] - 2*theta2[n]));

            while(ddtheta2[n]>=2*pi/(dt*dt))
                ddtheta2[n]=ddtheta2[n]-2*pi/(dt*dt);
            while(ddtheta2[n]<=-2*pi/(dt*dt))
                ddtheta2[n]=ddtheta2[n]+2*pi/(dt*dt);

            dtheta1[n]=dtheta1[n] + ddtheta1[n]*dt;
            while(dtheta1[n]>=2*pi/dt)
                dtheta1[n]=dtheta1[n]-2*pi/dt;
            while(dtheta1[n]<=-2*pi/dt)
                dtheta1[n]=dtheta1[n]+2*pi/dt;

            dtheta2[n]=dtheta2[n] + ddtheta2[n]*dt;
            while(dtheta2[n]>=2*pi/dt)
                dtheta2[n]=dtheta2[n]-2*pi/dt;
            while(dtheta2[n]<=-2*pi/dt)
                dtheta2[n]=dtheta2[n]+2*pi/dt;

            //cout<<n<<"\'"<<theta1[n]<<":"<<dtheta1[n]<<".."<<ddtheta1[n]<<"\n";

            if(dtheta1[n]>=2*pi/dt)
                exit(n);

            //drawing pendulum

            r = 255;
            gr = int(255-abs(dtheta2[n]/(2.*pi))*255);
            b = int(255-abs(dtheta1[n]/(2.*pi))*255);

            x1 = px0+2*floor(L1*sin(theta1[n]));
            y1 = py0+2*floor(L1*cos(theta1[n]));
            x2 = x1+2*floor(L2*sin(theta2[n]));
            y2 = y1+2*floor(L2*cos(theta2[n]));


            //double dt1 = dtheta1[n];
            //double dt2 = dtheta2[n];



            line(px0,py0,x1,y1,100,r,gr,b);
            line(x1,y1,x2,y2,100,r,gr,b);

        }
            string t = "test", dott = ".ppm", str;
            stringstream ss;
            ss<<i;
            ss>>str;
            t+=str; t+=dott;

            cout<<theta1[0]<<":"<<theta2[0]<<";;;"<<dtheta1[0]<<":"<<dtheta2[0]<<";;;"<<ddtheta1[0]<<":"<<ddtheta2[0]<<";;;"<<t<<":"<<"\n";

            saveppm(t);

    }
}
