#include "m_pd.h"

#define SQRT2 1.414213562373095048801688 

typedef struct _ihwt_tilde
{
    t_object x_obj;
    t_outlet *out1;
    t_float x_f;
} t_ihwt_tilde;

static t_class *ihwt_tilde_class;

static void *ihwt_tilde_new(t_floatarg t)
{
    t_ihwt_tilde *x = (t_ihwt_tilde *)pd_new(ihwt_tilde_class);
    outlet_new(&x->x_obj, &s_signal);
    x->x_f = 0;
    return (x);
}

static void ihwt_tilde_free(t_ihwt_tilde *x)
{
    outlet_free(x->x_obj.ob_outlet);
}

static t_int *ihwt_tilde_perform(t_int *w)
{
    t_ihwt_tilde *x = (t_ihwt_tilde *)(w[1]);
    t_sample  *in1 = (t_sample *)(w[2]);
    t_sample *out1 = (t_sample *)(w[3]);
    int          n = (int)(w[4]);

    t_sample *out1_ = out1;
    for (int i = 0; i < n; i++)
        *(out1_++) = *(in1++);

    int k = n;
    do
    {
        t_sample *p1 = out1, *p2 = out1 + (k >> 1);
        while (p2 <= out1_)
        {
            t_sample a = *p1, b = *p2;
            *p1 = (a + b)/SQRT2;
            *p2 = (a - b)/SQRT2;
            p1 += k;
            p2 += k;
        }
    } while ((k >>= 1) > 1);

    return w + 5;
}

static void ihwt_tilde_dsp(t_ihwt_tilde *x, t_signal **sp)
{
    dsp_add(ihwt_tilde_perform, 4, x,
        sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

void ihwt_tilde_setup(void)
{
    ihwt_tilde_class = class_new(gensym("ihwt~"),
        (t_newmethod)ihwt_tilde_new,
        (t_method)ihwt_tilde_free,
        sizeof(t_ihwt_tilde),
        CLASS_DEFAULT,
        A_DEFFLOAT, 0
    );
    class_addmethod(ihwt_tilde_class,
        (t_method)ihwt_tilde_dsp,
        gensym("dsp"),
        A_CANT, 0
    );
    CLASS_MAINSIGNALIN(ihwt_tilde_class, t_ihwt_tilde, x_f);
}