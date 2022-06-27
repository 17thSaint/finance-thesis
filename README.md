# Modeling a Financial System with Memory

In this project, I studied how memory can be included in the modeling of a financial system by using fractional calculus and fractional Brownian motion.  This research project was part of the Complex Systems profile during my MSc degree in Theoretical Physics at Utrecht University.  Brownian motion has long been used to model financial systems and fractional Brownian motion models have recently been used to exmaine the roughness of volatility.  However, these models all use analytical methods native to the financial sector.  This project used physics-based analysis tools including phase bevahior and dispersion relations, to study a financial model that includes memory.

## Fractional Calculus
There is no single method of extending regular calculus to fractional orders which means that there are a variety of definitions of the fractional derivative that are not always equivalent.  This project focuses on three of the most common versions: the Riemann-Liouville form, the Grünwald-Letnikov form, and the Caputo form.  The oldest of which is the Riemann-Liouville form whose derivative of order $\alpha$ is given by

$${}^{RL} D_{a+}^{\alpha} f(t) = \frac{1}{\Gamma(n-\alpha)} \frac{d^n}{dt^n} \int_{a}^{t} d\tau (t-\tau)^{n-\alpha-1} f(\tau)$$

which has a lower bound, a, and $\Gamma$ is the Gamma function: an extension of the factorial function.  The main issue, in this project, with the RL form is that as the integral reaches the upper bound, it begins to diverge numerically.  A discretized version of the RL fractional derivative is the Grünwald-Letnikov form given by

$${}^{GL} D_{a+}^{-\alpha} f(t) = \lim_{\substack{h\rightarrow0 \\ nh \rightarrow t-a}} h^{\alpha} \sum_{j=0}^{n}
\begin{pmatrix}
\alpha \\
j
\end{pmatrix}
f(t-jh)
$$

For this project, the Caputo form of the fractional derivative is used because it has the most effective inclusion of memory in its definition.  It is given by 

$${}^{C} D_{a+}^{\alpha} f(t) = \frac{1}{\Gamma(n-\alpha)} \int_{a}^{t} d\tau \frac{f^{(n)}(\tau)}{(t-\tau)^{\alpha-n+1}}, \quad n-1 < \alpha < n.$$

The Caputo form can also be calculated from the RL and GL forms using the following relation

$${}^{RL} D_{a}^{\alpha} f(t) = \sum_{k=0}^{n-1} \frac{f^{(k)}(a) (t-a)^{k-\alpha}}{\Gamma(k-\alpha+1)} + {}^{C} D_{a}^{\alpha} f(t).$$

## Fractional Brownian Motion

## The Financial Model

## Monte Carlo Solution Finder

## Results
