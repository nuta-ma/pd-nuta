#include "m_pd.h"
#include <math.h>
#include <stdint.h>
#define PI 3.14159265359
#define PI2 6.28318530718

typedef struct _morph_tilde
{
    t_object x_obj;
    t_inlet *in2;
    t_inlet *in3;
    t_inlet *in4;
    t_inlet *in5;
    t_inlet *in6;
    t_outlet *out1;
    t_outlet *out2;
    t_float x_ta;
    t_float x_tp;
    t_float x_f;
} t_morph_tilde;

static t_class *morph_tilde_class;

static void *morph_tilde_new(t_floatarg ta, t_floatarg tp)
{
    t_morph_tilde *x = (t_morph_tilde *)pd_new(morph_tilde_class);
    x->in2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->in3 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->in4 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    x->in5 = floatinlet_new(&x->x_obj, &x->x_ta);
    x->in6 = floatinlet_new(&x->x_obj, &x->x_tp);
    x->out1 = outlet_new(&x->x_obj, &s_signal);
    x->out2 = outlet_new(&x->x_obj, &s_signal);
    x->x_ta = ta;
    x->x_tp = tp;
    x->x_f = 0;
    return (x);
}

static void morph_tilde_free(t_morph_tilde *x)
{
    inlet_free(x->in2);
    inlet_free(x->in3);
    inlet_free(x->in4);
    inlet_free(x->in5);
    inlet_free(x->in6);
    outlet_free(x->out1);
    outlet_free(x->out2);
}

static double wrapphase(double c)
{
    c = fmod(c, PI2);
    return (PI <= c) ? (c - PI2) : (c < -PI) ? (c + PI2) : c;
}

static double avgphase(double a, double b, double t)
{
    return a + t * wrapphase(b - a);
}


static t_int *morph_tilde_perform(t_int *w)
{
    t_morph_tilde *x = (t_morph_tilde *)(w[1]);
    t_sample *aa = (t_sample *)(w[2]);
    t_sample *pa = (t_sample *)(w[3]);
    t_sample *ab = (t_sample *)(w[4]);
    t_sample *pb = (t_sample *)(w[5]);
    t_sample *ac = (t_sample *)(w[6]);
    t_sample *pc = (t_sample *)(w[7]);
    int N = (int)(w[8]), n = N/2 + 1;

    t_sample ta = (x->x_ta < 0) ? 0 : (1 < x->x_ta) ? 1 : x->x_ta;
    t_sample tp = (x->x_tp < 0) ? 0 : (1 < x->x_tp) ? 1 : x->x_tp;

    double *iaa = (double *)getbytes((n + 1) * sizeof(double));
    double *iab = (double *)getbytes((n + 1) * sizeof(double));
    double *iac = (double *)getbytes(n * sizeof(double));
    double *pc1 = (double *)getbytes(n * sizeof(double));
    double *pc2 = (double *)getbytes(n * sizeof(double));

    unsigned int i = 0, j = 0, k, l;

    iaa[0] = -aa[0] / 2;
    iab[0] = -ab[0] / 2;

    for (k = 0; k < n; k++)
    {
        iaa[k + 1] = iaa[k] + aa[k];
        iab[k + 1] = iab[k] + ab[k];
    }

    iaa[n] -= aa[n - 1] / 2;
    iab[n] -= ab[n - 1] / 2;

    double scale = (1 - ta) * iaa[n] + ta * iab[n];

    for (k = 0; k < n; k++)
    {
        iaa[k] /= iaa[n];
        iab[k] /= iab[n];
    }

    iaa[n] = 1 + aa[n - 1] / 2 / iaa[n];
    iab[n] = 1 + ab[n - 1] / 2 / iab[n];

    double t, f, ia, p1, p2, f_prev = 0, ia_prev = 0, p1_prev = 0, p2_prev = 0;
    do
    {
        if (iab[j + 1] < iaa[i + 1])
        {
            j++;
            t = (iab[j] - iaa[i]) / (iaa[i + 1] - iaa[i]);
            f = pow((i + t) / j, 1 - ta) * j;
            ia = iab[j];
            p1 = avgphase((i == 0) ? avgphase(-pa[1], pa[0], t) : avgphase(pa[i - 1], pa[i], t), pb[j - 1], tp);
            p2 = avgphase((i == n - 1) ? avgphase(pa[n - 1], -pa[n - 2], t) : avgphase(pa[i], pa[i + 1], t), (j == n) ? -pb[n - 2] : pb[j], tp);
        }
        else
        {
            i++;
            t = (iab[j] == iab[j + 1]) ? 0 : (iaa[i] - iab[j]) / (iab[j + 1] - iab[j]);
            f = i * pow((j + t) / i, ta);
            ia = iaa[i];
            p1 = avgphase(pa[i - 1], (j == 0) ? avgphase(-pb[1], pb[0], t) : avgphase(pb[j - 1], pb[j], t), tp);
            p2 = avgphase((i == n) ? -pa[n - 2] : pa[i], (j == n - 1) ? avgphase(pb[n - 1], -pb[n - 2], t) : avgphase(pb[j], pb[j + 1], t), tp);
        }

        for (k = (unsigned int) f_prev; k < (unsigned int) f; k++)
        {
            t = (k + 1 - f_prev) / (f - f_prev);
            iac[k] = ia_prev + (ia - ia_prev) * t;
            pc1[k] = avgphase(p1_prev, p1, t);
            pc2[k] = avgphase(p2_prev, p2, t);
        }

        f_prev = f;
        ia_prev = ia;
        p1_prev = p1;
        p2_prev = p2;
    }
    while (ia < 1);

    ac[0] = 2 * scale * iac[0];
    pc[0] = wrapphase(pc1[0]);

    for (l = 1; l < k; l++)
    {
        ac[l] = scale * (iac[l] - iac[l - 1]);
        pc[l] = wrapphase(avgphase(pc1[l], pc2[l - 1], 0.5));
    }

    if (l < n - 1)
    {
        ac[l] = scale * (1 - iac[l - 1]);
        pc[l] = wrapphase(avgphase(pc1[l], pc2[l - 1], 0.5));
    }
    else
    {
        ac[n - 1] = 2 * scale * (1 - iac[n - 2]);
        pc[n - 1] = wrapphase(pc2[n - 2]);
    }

    l++;

    while (l < N)
    {
        ac[l] = 0;
        pc[l] = 0;
        l++;
    }

    freebytes(iaa, (n + 1) * sizeof(double));
    freebytes(iab, (n + 1) * sizeof(double));
    freebytes(iac, n * sizeof(double));
    freebytes(pc1, n * sizeof(double));
    freebytes(pc2, n * sizeof(double));

    return (w + 9);
}

static void morph_tilde_dsp(t_morph_tilde *x, t_signal **sp)
{
    dsp_add(morph_tilde_perform, 8, x,
        sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec,
        (t_int)sp[0]->s_n);
}

void morph_tilde_setup(void)
{
    morph_tilde_class = class_new(gensym("morph~"),
        (t_newmethod)morph_tilde_new,
        (t_method)morph_tilde_free,
        sizeof(t_morph_tilde),
        CLASS_DEFAULT,
        A_DEFFLOAT,
        A_DEFFLOAT, 0
    );
    class_addmethod(morph_tilde_class,
        (t_method)morph_tilde_dsp,
        gensym("dsp"),
        A_CANT, 0
    );
    CLASS_MAINSIGNALIN(morph_tilde_class, t_morph_tilde, x_f);
}