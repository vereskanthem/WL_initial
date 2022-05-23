#include <math.h>
#include <stdio.h>

int T_step_max = 800;
double T0 = 0.1, T1 = 8.1;

int main()
{
    int Lx,Ly,Lz, N, n;
    int i, T_step;
    double prom, *ln_g, *H, *E;
    double E_eq, T, dT, Z, E2_eq;
    char fname[80];
    FILE *OUT;

    OUT = fopen("temp.cfg","r");
    fscanf(OUT, "%d", &Lx);
    fscanf(OUT, "%d", &Ly);
    fscanf(OUT, "%d", &Lz);
    fscanf(OUT, "%d", &n);
    fclose(OUT);

    N = Lx*Ly*Lz;

    ln_g = new double [n];
    H = new double [n];
    E = new double [n];

    dT = (T1-T0)/T_step_max;

    OUT = fopen("data.dat","r");

    for(i=0;i<n;++i)
    {
        fscanf(OUT,"%lf", &E[i]);
        fscanf(OUT,"%lf", &ln_g[i]);
        fscanf(OUT,"%lf", &H[i]);
    }
    fclose(OUT);

    OUT = fopen("result.dat","w+");
    for(T_step=0;T_step<=T_step_max;++T_step)
    {
        T = T0 + dT*T_step;
        fprintf(OUT,"%.4lf", T);

        prom = ln_g[0] - E[0]/T;

        for(i=1;i<n;++i)
        {
            if((ln_g[i] - E[i]/T)>prom)
            {
                prom = ln_g[i] - E[i]/T;
            }
        }

        if(T_step%10==0)
        {
            printf("T = %.5lf", T);
            printf("\n");
        }

        E_eq = 0.0;
        E2_eq = 0.0;
        Z = 0.0;

        for(i=0;i<n;++i)
        {
            E_eq = E_eq + E[i]*exp(ln_g[i]-E[i]/T - prom);
            E2_eq = E2_eq + E[i]*E[i]*exp(ln_g[i]-E[i]/T - prom);
            Z = Z + exp(ln_g[i]-E[i]/T - prom);
        }
        fprintf(OUT, "\t%.9lf", (E_eq/Z)/N);
        fprintf(OUT, "\t%.9lf", (((E2_eq/Z)-(E_eq/Z)*(E_eq/Z))/(T*T))/N);
        fprintf(OUT, "\t%.9lf", (-T*prom-(T)*log(Z))/N);
        fprintf(OUT, "\t%.9lf", (((E_eq/Z)-(-T*prom-(T)*log(Z)))/T)/N);
        fprintf(OUT, "\n");
    }
    fclose(OUT);

    delete(H);
    delete(ln_g);
    delete(E);
    return 0;
}
