import numpy as np
import pandas as pd
import scipy.optimize
import scipy.stats as st
import scipy.special
import numba



def log_like_gamma(params, n):
    """Log likelihood for gamma measurements."""
    α, β = params
    
    if α <= 0 or β <= 0:
        return -np.inf

    return np.sum(st.gamma.logpdf(n, α, scale=1/β))



def mle_gamma(n):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        res = scipy.optimize.minimize(
            fun=lambda params, n: -log_like_gamma(params, n),
            x0=np.array([1, 1]),
            args=(n,),
            method="Powell",
        )

    if res.success:
        α_mle, β_mle = res.x
    else:
        raise RuntimeError("Convergence failed with message", res.message)

    return α_mle, β_mle



def log_like_custom_model(params, n):
    β1, β2 = params
    
    #Set physical constraints
    if β1 <= 0 or β2 <= 0:
        return -np.inf
    
    #Constrain further to β1 < β2 since exp[β1 * t] > exp[β2 * t]
    if β1 > β2:
        return -np.inf
    
    if np.isclose(β1, β2):
        return np.sum(-n * β2 + np.log(n) + np.log(β2**2)) #from limit calculator
    
    #otherwise, use derived log-likelihood
    log_model = (np.log((β1 * β2) / (β2 - β1))
                 - β1 * n 
                 + np.log(1 - np.e**-((β2 - β1) * n))
                )
    
    return np.sum(log_model)



def mle_model(n):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        res = scipy.optimize.minimize(
            fun=lambda params, n: -log_like_custom_model(params, n),
            x0=np.array([1, 2]),
            args=(n,),
            method="Powell",
        )

    if res.success:
        β1_mle, β2_mle = res.x
    else:
        raise RuntimeError("Convergence failed with message", res.message)

    return β1_mle, β2_mle




def custom_mle(n):
    """Define our own mle"""
    return np.array([2 / n.mean()]*2)




def gamma_sample(x, α, β):
    """CDF for gamma distribution"""
    return st.gamma.cdf(x, α, scale=1 / β)




def model_sample(x, β1, β2):
    """Sample the betas with the model"""
    if β1 == β2:
        return -np.inf
    else:
        return (β1 * β2 / (β2 - β1)) * (
            1 / β1 * (1 - np.exp(-β1 * x)) - 1 / β2 * (1 - np.exp(-β2 * x))
        )



def comp_plot(df):
    """Compare the data to the model as ECDFs"""
    n = df["time to catastrophe (s)"].values

    compare = iqplot.ecdf(
        df, q="time to catastrophe (s)", style="staircase", conf_int=True
    )

    x = np.linspace(n.min(), n.max(), 400)
    y_gamma = gamma_sample(x, *mle_gamma(n))
    y_model = model_sample(x, *mle_model(n))

    compare.line(x, y_gamma, color=colorcet.b_glasbey_category10[1], width=2, alpha=0.8)
    compare.line(x, y_model, color=colorcet.b_glasbey_category10[2], width=2, alpha=0.8)

    legend = bokeh.models.Legend(
        items=[
            (
                "experiment",
                [compare.line(color=colorcet.b_glasbey_category10[0], width=2)],
            ),
            (
                "gamma",
                [
                    compare.line(
                        color=colorcet.b_glasbey_category10[1], width=2, alpha=0.8
                    )
                ],
            ),
            (
                "5.2/8.2 model",
                [
                    compare.line(
                        color=colorcet.b_glasbey_category10[2], width=2, alpha=0.8
                    )
                ],
            ),
        ],
        location="center",
    )

    compare.add_layout(legend, "right")

    bokeh.io.show(compare)

    return compare



def aic(df):

    n = df["time to catastrophe (s)"].values

    AIC_gamma = -2 * log_like_gamma(mle_gamma(n), n)
    AIC_model = -2 * log_like_model(mle_model(n), n)

    w_gamma = (np.exp(-(AIC_gamma - max(AIC_gamma, AIC_model)) / 2)) / (
        np.exp(-(AIC_gamma - max(AIC_gamma, AIC_model)) / 2)
        + np.exp(-(AIC_model - max(AIC_gamma, AIC_model)) / 2)
    )

    w_model = (np.exp(-(AIC_model - max(AIC_gamma, AIC_model)) / 2)) / (
        np.exp(-(AIC_gamma - max(AIC_gamma, AIC_model)) / 2)
        + np.exp(-(AIC_model - max(AIC_gamma, AIC_model)) / 2)
    )

    print("Akaike weight for gamma = " + str(w_gamma))
    print("Akaike weight for 5.2/8.2 model = " + str(w_model))
    print("Gamma distribution is favored by a factor of " + str(w_gamma / w_model))

    return w_gamma, w_model, w_gamma / w_model



def get_conc(df):
    """Get unique concentration values for the data frame"""
    return [val for val in df["concentration (μM)"].unique()]



def get_param(conc):
    """Return and print the parameter values"""
    α = mle_gamma(df.loc[df["concentration (μM)"] == conc])[0]
    β = mle_gamma(df.loc[df["concentration (μM)"] == conc])[1]

    print(
        str(conc)
        + "μM:"
        + (2 - len(str(conc))) * " "
        + " α = "
        + str(mle_gamma(df.loc[df["concentration (μM)"] == conc])[0])
        + "\n      β = "
        + str(mle_gamma(df.loc[df["concentration (μM)"] == conc])[1])
    )

    return α, β




