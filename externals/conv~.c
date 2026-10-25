#include "m_pd.h"
#include <math.h>
#include <string.h>
#include <pthread.h>

#define PI 3.1415926535897932384626433

typedef struct conv_arg_{
    t_sample *c_re; /* carrier */
    t_sample *c_im;
    t_sample *i_re; /* impulse response */
    t_sample *i_im;
    int N; /* 2n */
} conv_arg_t;

/* unsorted fast Fourier transform */
static void ufft(t_sample *re, t_sample *im, int N)
{
    int M = 1, s;
    t_sample lr, li, rr, ri, a, b, c, d;
    while (N!=1)
    {
        N /= 2;
        for (int i=0;i<M;i++) for(int j=0;j<N;j++)
        {
            s = 2*i*N+j;
            lr = re[s];
            li = im[s];
            rr = re[s+N];
            ri = im[s+N];

            re[s] = lr+rr;
            im[s] = li+ri;

            a = lr-rr;
            b = li-ri;
            c = cos(PI*j/N);
            d = sin(PI*j/N);

            re[s+N] =  a*c + b*d;
            im[s+N] = -a*d + b*c;
        }
        M *= 2;
    }
}

/* inverse unsorted fast Fourier transform */
static void iufft(t_sample *re, t_sample *im, int M)
{
    int N = 1, s;
    t_sample a, b, c, d, lr, li, rr, ri;
    while (M!=1)
    {
        M /= 2;
        for (int i=0;i<M;i++) for(int j=0;j<N;j++)
        {
            s = 2*i*N+j;
            a = re[s+N];
            b = im[s+N];
            c = cos(PI*j/N);
            d = sin(PI*j/N);

            lr = re[s];
            li = im[s];
            rr = a*c - b*d;
            ri = a*d + b*c;

            re[s] = lr+rr;
            im[s] = li+ri;

            re[s+N] = lr-rr;
            im[s+N] = li-ri;
        }
        N *= 2;
    }
    for(int i=0;i<N;i++)
    {
        re[i] /= N;
        im[i] /= N;
    }
}

/* multiply complex numbers (Hadamard product) */
static void multi(t_sample *re, t_sample *im, t_sample *re_, t_sample *im_, int N)
{
    for(int i=0;i<N;i++){
        re[i] = re[i]*re_[i] - im[i]*im_[i];
        im[i] = re[i]*im_[i] + im[i]*re_[i];
    }
}

/* convolution */
static void conv(conv_arg_t *a)
{
/*    ufft(a->c_re, a->c_im, a->N);
    ufft(a->i_re, a->i_im, a->N);
    multi(a->c_re, a->c_im, a->i_re, a->i_im, a->N);
    iufft(a->c_re, a->c_im, a->N); */
}

static t_class *conv_tilde_class;

typedef struct _conv_tilde{
    t_object x_obj;
    t_float f;
    t_inlet *in2;
    t_sample *pr; /* previous block */
    int npr;
    int ph; /* phase: which pair to convolute */
    pthread_t th[2];
    conv_arg_t a[2];
} t_conv_tilde;

static t_int *conv_tilde_perform(t_int *w)
{
    t_conv_tilde *x = (t_conv_tilde *)(w[1]);
    t_sample   *in1 =     (t_sample *)(w[2]);
    t_sample   *in2 =     (t_sample *)(w[3]);
    t_sample   *out =     (t_sample *)(w[4]);
    int           n =            (int)(w[5]);
    conv_arg_t   *a =           &x->a[x->ph];

    if(pthread_join(x->th[x->ph], NULL))
    {
        post("joinできない");
    }
    else
    {
        int i = n;
        while(i--) out[n-1-i] = in1[n-1-i]*0.5;//a->c_re[a->N-1-i];

    }

    freebytes(a->c_re, a->N*sizeof(t_sample));
    freebytes(a->c_im, a->N*sizeof(t_sample));
    freebytes(a->i_re, a->N*sizeof(t_sample));
    freebytes(a->i_im, a->N*sizeof(t_sample));
    a->N    = 2*n;
    a->c_re = getbytes(a->N*sizeof(t_sample));
    a->c_im = getbytes(a->N*sizeof(t_sample));
    a->i_re = getbytes(a->N*sizeof(t_sample));
    a->i_im = getbytes(a->N*sizeof(t_sample));

    if(x->npr>n)
    {
        x->npr = n;
        x->pr += x->npr - n;
    }

    
    memset(         a->c_re,     0, (n-x->npr)*sizeof(t_sample));
    memcpy(n-x->npr+a->c_re, x->pr,     x->npr*sizeof(t_sample));
    memcpy(       n+a->c_re,   in1,          n*sizeof(t_sample));
    memset(         a->c_im,     0,       a->N*sizeof(t_sample));
    memcpy(         a->i_re,   in2,          n*sizeof(t_sample));
    memset(       n+a->i_re,     0,          n*sizeof(t_sample));
    memset(         a->i_im,     0,       a->N*sizeof(t_sample));

    if(pthread_create(&x->th[x->ph], NULL, (void *)conv, (void *)a)) post("createできない");

    x->ph  = !x->ph;
    x->pr  = (t_sample*)resizebytes(x->pr, n, n*sizeof(t_sample));
    x->npr = n;

    return w + 6;
}

static void conv_tilde_dsp(t_conv_tilde *x, t_signal **sp)
{
    dsp_add(conv_tilde_perform, 5, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n);
}

static void dummy(char *a)
{
}

static void *conv_tilde_new(void)
{
    post("よお");
    t_conv_tilde *x = (t_conv_tilde *)pd_new(conv_tilde_class);

    x->f = 0;

    x->pr  = (t_sample *)getbytes(sizeof(t_sample));
    x->npr = 1;

    x->ph = 0;
    conv_arg_t *a;
    for(int i = 0; i<2; i++)
    {
        a = &x->a[i];
        a->c_re = getbytes(4*sizeof(t_sample));
        a->c_im = getbytes(4*sizeof(t_sample));
        a->i_re = getbytes(4*sizeof(t_sample));
        a->i_im = getbytes(4*sizeof(t_sample));
        a->N = 4;
       if(pthread_create(&x->th[i],NULL,(void *)dummy,(void *)"a")) post("createできない");
    }
    x->in2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    outlet_new(&x->x_obj, &s_signal);
    return (void *)x;
}

static void conv_tilde_free(t_conv_tilde *x)
{
    inlet_free(x->in2);
    outlet_free(x->x_obj.ob_outlet);

    freebytes(x->pr, x->npr*sizeof(t_sample));

    conv_arg_t *a;
    for(int i = 0; i<2; i++)
    {
        a = &x->a[i];
        freebytes(a->c_re, a->N*sizeof(t_sample));
        freebytes(a->c_im, a->N*sizeof(t_sample));
        freebytes(a->i_re, a->N*sizeof(t_sample));
        freebytes(a->i_im, a->N*sizeof(t_sample));

        pthread_detach(x->th[i]);
    }
}

void conv_tilde_setup(void)
{
    conv_tilde_class = class_new(
        gensym("conv~"),
        (t_newmethod)conv_tilde_new,
        (t_method)conv_tilde_free,
        sizeof(t_conv_tilde),
        CLASS_DEFAULT,
        0);
    class_addmethod(conv_tilde_class,
        (t_method)conv_tilde_dsp,
        gensym("dsp"),
        A_CANT, 0);
    CLASS_MAINSIGNALIN(conv_tilde_class, t_conv_tilde, f);
}