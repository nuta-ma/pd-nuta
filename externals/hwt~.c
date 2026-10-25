#include "m_pd.h"

#define SQRT2 1.414213562373095048801688 

typedef struct _hwt_tilde
{
    t_object x_obj;
    t_outlet *out1;
    t_float x_f;
} t_hwt_tilde;

static t_class *hwt_tilde_class;

static void *hwt_tilde_new(t_floatarg t)
{
    t_hwt_tilde *x = (t_hwt_tilde *)pd_new(hwt_tilde_class);
    outlet_new(&x->x_obj, &s_signal);
    x->x_f = 0;
    return (x);
}

static void hwt_tilde_free(t_hwt_tilde *x)
{
    outlet_free(x->x_obj.ob_outlet);
}

static t_int *hwt_tilde_perform(t_int *w)
{
    t_hwt_tilde *x = (t_hwt_tilde *)(w[1]);
    int          n = (int)(w[4]);
    t_sample  *in1 = (t_sample *)(w[2]) + n - 1;
    t_sample *out1 = (t_sample *)(w[3]) + n - 1;

    for (int i = 1; i < n; i += 2)
    {
        *(out1 - i + 1) = *(in1 - i + 1);
        *(out1 - i) = *(in1 - i);

        int k = 1;
        do
        {
            t_sample a = *(out1 - i), b = *(out1 - (i&~k));
            *(out1 - i)      = (a + b)/SQRT2;
            *(out1 - (i&~k)) = (a - b)/SQRT2;
        } while ((i & (k <<= 1)) != 0);
    }

    return w + 5;
}

static void hwt_tilde_dsp(t_hwt_tilde *x, t_signal **sp)
{
    dsp_add(hwt_tilde_perform, 4, x,
        sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

void hwt_tilde_setup(void)
{
    hwt_tilde_class = class_new(gensym("hwt~"),
        (t_newmethod)hwt_tilde_new,
        (t_method)hwt_tilde_free,
        sizeof(t_hwt_tilde),
        CLASS_DEFAULT,
        A_DEFFLOAT, 0
    );
    class_addmethod(hwt_tilde_class,
        (t_method)hwt_tilde_dsp,
        gensym("dsp"),
        A_CANT, 0
    );
    CLASS_MAINSIGNALIN(hwt_tilde_class, t_hwt_tilde, x_f);
}