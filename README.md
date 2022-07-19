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
Fractional Brownian motion can be generated using two different methods.  The first method uses a differential equation called the Langevin equation and comes from a physics-based method.  The second comes from the financial world and uses a stochastic process to model the motion.

### Langevin Equation
In order to describe the physics approach to fractional Brownian motion, we derive the Langevin equation:

$$
M \ddot{x} + \gamma \dot{x} = f(t),
$$

where $f(t)$ is regular Gaussian noise.  Gaussian noise is the time deriative of regular Browian motion.  The resulting motion $x(t)$ is a Browian motion path of our trace particle.  This derivation is done by starting with the Caldeira-Leggett model of quantum dissipation.  The derivation starts by building a Lagrangian for the model which comes in four parts: the trace particle of our system, $\mathcal{L_S}$, the bath, $\mathcal{L_B}$, the interaction, $\mathcal{L_I}$, and a counter-term, $\mathcal{L_{CT}}$.  The full details are given in the report.  In order to get fractional Brownian motion, the Lagrangian can be solved to find the equations of motion for a spectral function $J(\omega)$.  When the spectral function goes as a power law given by:

$$
J(\omega) = \gamma \frac{\pi}{2} \frac{\sec(\pi \alpha / 2)}{\Gamma(\alpha) \Gamma(1-\alpha)} \omega^{\alpha},
$$

the resulting fractional Langevin equation now has the friction term as a fractional derivative of order $\alpha$ instead of the first-order derivative seen in the regular Langevin equation.

$$
M \frac{d^2x}{dt^2} + \gamma {}^{C}D_{0}^{\alpha} \big[ x(t) \big] = f(t).
$$

The time correlation of the noise for this fractional Langevin equation is thus given by

$$
\big\langle f(t) f(t') \big\rangle = \frac{\gamma k_B T}{\Gamma(1-\alpha)} \vert t-t' \vert^{-\alpha}. 
$$

### Stochastic Process
In stochastic language, regular Brownian motion, $B(t)$, is a continuous-time process with independent increments: a Wiener process.  A process of this kind has increments that are randomly distributed and scaled according to the size of the increments.  The number of increments in the process remains constant for depending on the time-scale of the process the size of the increments must be rescaled.  This random distribution is what leads the time correlation to be a Dirac delta.  However, in stochastics language, the time correlation is called the time covariance because it is the covariance of a motion at two points in time; from now on it will be referred to simply as the covariance.

In order to include memory in the system, these increments are changed such that successive increments are dependent on those that came before.  The parameter that defines the strength of that historical dependence is called the Hurst parameter, $H$, where $H$ = 1/2 is regular Brownian motion.  With this, the covariance of a fractional Brownian motion path, $B_H (t)$, is given by

$$
\big\langle B_{H}(t) B_{H}(\tau) \big\rangle = \frac{\sigma_{0}^2 \Gamma(2-2H)}{4H\Gamma(3/2-H)\Gamma(1/2+H)} \bigg( t^{2H} + \tau^{2H} - \vert t-\tau \vert^{2H} \bigg),
$$

where $\sigma_0$ is a constant that comes from regular Brownian motion.  From here on, this value is normalized out of the paths.  The Hurst parameter is limited to $H \in (0,1)$, and the question is what these different values of H actually correspond to.

include history correlation figure

Figure here shows the time autocorrelation of fractional Brownian motion paths as a function of $H$.  The time autocorrelation measures the correlation between a point in the path at time $t$ with a point on the path some time later.  The correlation is averaged over points from the entire given path.  What is seen, in Figure here, is that for $H$ = $1/2$ the correlation to history is zero, i.e. no memory.  For $H$ $<$ $1/2$, there is negative correlation to history while for $H$ $>$ $1/2$, there is a positive correlation.  We know that this correlation as a function of the Hurst parameter is exactly correct.  Also, as expected, the strength of the correlation decays as the time separation gets larger.

However, it is important to know that motion is not the same as noise.  The increment process $\xi_H (t)=B_H (t+1) - B_H (t)$ is fractional Gaussian noise.  Now, because we are looking at a continuous-time process this expression needs to be continuous in time and so we get $\xi_H (t) = d B_H / dt$.  Above it was mentioned that $H$ = 1/2 is regular Brownian motion, so the noise generated by this motion is white noise.  The expression given in Equation eq describes the covariance of the motion, not the noise.  The noise covariance is given by: 

$$
\langle \xi_{H}(t) \xi_{H}(\tau) \rangle = \frac{ (2H-1)\Gamma(2-2H)}{2\Gamma(3/2-H)\Gamma(1/2+H)} \vert t - \tau \vert^{2H-2}.
$$

The generalization of white noise generates what is known as coloured noise.  This power law relation also looks quite similar to the time correlation from the Langevin equation approach.

### Equivalence of Approaches
The aim of the previous sections was to describe the two different approaches to generating fractional Brownian motion.  Theoretically, they both generate the same dynamics; hence they should be related to each other.  There is an obvious similarity in the power law form between time correlation of the two approaches.  From the physics point of view, it is known that $f(t)$ is some kind of noise and from the financial stochastic point of view, it is known that $\xi_H (t)$ is noise.  Let's set them equal to each other, up to a constant $\epsilon$.  This means that we have the relation

$$
\boxed{\alpha = 2-2H,}
$$

Another point needs to be addressed as well: the phrase "up to a constant".  There are some coefficients in front of the time correlation equations that are non-trivial.  Taking the noise $f(t)$ to be the colored noise generated from a stochastic process means that the fractional Langevin equation from becomes

$$
M \frac{d^2x}{dt^2} + \gamma {}^{C}D_{0}^{\alpha} \big[ x(t) \big] = \epsilon \xi_H (t)
$$

along with the caveat of the coefficients:

$$
\epsilon^2 = 2\gamma k_B T \frac{\Gamma(3/2-H)\Gamma(1/2+H)}{\Gamma(2H)\Gamma(2-2H)}.
$$

The word dissipation was used to describe the transfer of energy between the system and the bath and this was an intentional choice of words.  The coefficient $\epsilon$, maintains the fluctuation-dissipation theorem discussed in the Langevin case and connects to the stochastics approach as well.

## The Financial Model
### Analytic Description

Here, we give a bird's eye view of how a financial market can be modeled using the Langevin equation based on a model by Cont and Bouchaud.  It begins with the dynamics of simple supply and demand.  Assuming that supply and demand start off in equilibrium, for constant supply, when demand increases, the price increases because of the competition between actors.  Then, for constant demand, when supply increases, the price decreases because the value is tied to the availability.  Therefore, there exists a direct relationship between price and demand and an inverse relationship between price and supply.  For price, $x(t)$, demand $\phi_{D}$, and supply $\phi_{S}$,

$$
\frac{dx}{dt} = u(t) = \frac{\Delta \phi}{\lambda}, \quad \Delta \phi = \phi_{D} - \phi_{S},
$$

where $u(t)$ are the returns and $\lambda$ is called the market depth.  The market depth is the amount of excess demand required to move the price by one unit, or in financial terms: market liquidity, from now on referred to just as liquidity.

figure of liquidity table

In a financial market, there is an orderbook: a list of all prices that the actors are willing to buy and sell the asset.  When that orderbook only has a small number of actors willing to buy/sell the asset at only a few prices, as in the right side of the table, the market is illiquid.  When large transactions are made, the orderbook will empty out and the current price of the asset will change drastically.  When the orderbook is well filled, as in the left side of the table, there are large numbers of actors willing to buy and sell at any price around the current value: high liquidity.  This means that even when a transaction is made that buys/sells a large amount of the asset, the orderbook stays well filled and mostly unchanged.

Intuitively, the market depth, $\lambda$, is analogous to mass.  Large masses have a lot of inertia and thus are less subject to changes in their position.  This is mirrored in $\lambda$, where substantial market depth means that the price is less subject to change.

Market makers are the entities that satisfy orders for actors, and the rate at which they can fulfill those orders is given by the parameter $\Gamma$, which is large for highly liquid markets.  We assume that market makers act symmetrically such that they complete buy orders at the same rate as sell orders.  This allows us to expand $\Gamma$ as

$$
\Gamma (\phi) = \beta + \beta' \phi + ...
$$

where $\beta$ is a lowest-order approximation to the rate at which market makers satisfy orders.  The equation for the market makers is then:

$$
\frac{d\phi_{D,S}}{dt} \bigg\vert_{MM} = -\Gamma(\phi_{S,D}) \phi_{D,S},
$$
$$
\frac{d \Delta \phi}{dt} \bigg\vert_{MM} = \lambda \frac{d^2 x}{dt^2} = -\beta \big( \phi_D - \phi_S \big) = - \beta \lambda \frac{dx}{dt}.
$$

Now we have a differential equation that describes the change in price given that the market makers process transactions at the rate $\Gamma$.

The change in price of an asset, from the actor point of view, can be modeled as a function of two terms: a random term, $f_{D,S}$(t), and a trend term, $m_{D,S}$(t).

$$
\frac{d \phi_{D,S}}{dt} \bigg\vert_A = m_{D,S} (t) + f_{D,S} (t).
$$

The trend term is a function of the anticipated return, R(t), and the anticipated risk, $\Sigma$(t), for demand and supply.  In this particular case, we assume that the asset is risk neutral, which is a first-order approximation in the returns that assumes the risk to be constant in time: $\Sigma$(t)=$\Sigma^0$.  The anticipated return is given by

$$
R_{D,S}(t) = R_{D,S}^{0} + a_{D,S} \int_{0}^{t} d \tau K_R (t-\tau) \frac{dx}{d \tau},
$$

where $a_{D,S}$ measures the impact of the observed recent trend on the anticipated return and $K_R$ is a normalized kernel that describes how the historical trend is measured by the users.  We assume that the asset price has no general trend up or down, so the price is only determined by the demand/supply interaction and the random term.  This mean-reverting behavior is a reasonable assumption because the system here is isolated and thus is not subject to external factors.    It then follows that the constant terms in the risk and return cancel: $R_{D,S}^{0}$=$\Sigma_{D,S}^{0}$.

The form of the anticipated return looks very similar to a fractional derivative.  Indeed if we take

$$
K_R(t)=\frac{t^{-\alpha}}{\Gamma(1-\alpha)},
$$

we get the Caputo fractional derivative of price.  Putting together the return and risk for demand and supply, we get:

$$
\frac{d \Delta \phi}{dt} \bigg\vert_{A} = \frac{(a_D - a_S)}{\Gamma(1-\alpha)} \int_{0}^{t} d \tau (t-\tau)^{-\alpha} \frac{dx}{dt} + \big( R_{D}^{0} - \Sigma_{D}^{0} \big) - \big( R_{S}^{0} - \Sigma_{S}^{0} \big) + f(t), \quad \textrm{and}
$$
$$
\lambda \frac{d^2 x}{dt^2} \bigg\vert_A = a D^{\alpha} \big[ x(t) \big] + f(t),
$$

where we use the null-trend condition discussed above to move to the second line, $f(t)$ = $f_D$(t) - $f_S$(t), and $a$ = $a_D - a_S$.

To get the full stochastic differential equation for price, we combine the market maker component and the actor component to get:

$$
\lambda \frac{d^2 x}{dt^2} + \beta \lambda \frac{dx}{dt} - a D^{\alpha} \big[ x(t) \big]  = f(t).
$$

In the limiting case of no memory, we know that the fractional derivative term in the differential equation just goes to a first order derivative like a regular friction term.  Therefore, we take $\alpha$ = 1 as an upper bound and allow $\alpha$ to range in $(0,1)$.  For $\alpha$ in this range, that means that the Hurst parameter, $H \in (1/2,1)$.  As discussed in the beginning, when $H$ $>$ $1/2$, the path has a positive correlation with history.  From a financial point of view this actually makes a lot of sense.  If $H$ $<$ $1/2$, then there would be a negative correlation with history meaning that when the price goes up, actors expect the price to go down next.  This idea is quite contrary to the fundamental behavior of financial systems so it is reasonable to restrict $\alpha$ to this range.

### Noise and the Fluctuation-Dissipation Theorem

We know the covariance of the random term as a function of the market depth, $\lambda$, and a parameter $\Theta$:

$$
\langle f(t) f(t') \rangle = 2 \lambda^2 \Theta \delta (t-t')
$$

where $\Theta$ is the susceptibility of the market to random fluctuations.  The whole aim of our description here is to look at the impact of history, and in principle, $\Theta$ should be dependent on history.  If the system was hit by a large event recently, actors are on edge and will be more reactionary to incoming external news.  So, the covariance will be,

$$
\langle f(t) f(t') \rangle = 2 \lambda^2 \big[ \Theta_1 \delta (t-t') + \Theta_2 (t-t')^{-\alpha} \big].
$$

Therefore, the random term, $f(t)$, can be thought of as the sum of two, uncorrelated, noise terms.  One of them is white noise, which leads to the $\delta(t-t')$, and the other is a colored noise term which leads to the $(t-t')^{-\alpha}$ because of the memory.

Mathematically, this also makes sense when we look back at differential equation, which has an integer-order derivative and a fractional-order derivative.  Using the fluctuation-dissipation theorem, the integer-order term interacts with the white noise and the fractional-order term interacts with the colored noise.  Then, we notice that the fractional-order term is negative. Thus, the random term $f(t)$ should be written as

$$
f(t) = \eta(t) - \xi_{\alpha} (t),
$$

where $\eta$ is the white noise and $\xi_{\alpha}$ is the colored noise.  Both of the noises include constants that satisfy the fluctuation-dissipation theorem: $b_W$ and $b_C$.  Using all of these relationships, we find the following relations between the constants,

$$
\Theta_1 = \frac{\beta k_B T}{\lambda}, \quad b_{W}^{2} = 2 \beta \lambda k_B T; \qquad \Theta_2 = \frac{a k_B T}{2 \lambda^2 \Gamma(1-\alpha)}, \quad b_{C}^{2} = a k_B T \frac{2 \Gamma((1+\alpha)/2) \Gamma((3-\alpha)/2)}{\Gamma(2-\alpha)\Gamma(\alpha)}.
$$

These coefficients all include a factor of $k_B T$ which in physics language is the thermal energy of the system.  In physics, thermal energy corresponds to the average velocity of particles within a system.  In a financial system, position can be thought of as the amount of an asset an actor has and thus velocity is how much the actor buys/sells the asset.  The thermal energy then could be related to the trading volume at a given time.  High temperature financial systems would be ones that are traded frequently while low temperature ones correspond to quieter systems.

With these two components, we essentially have a particle, our price, interacting with two different, uncorrelated, baths.  The first bath has no memory and is just pushing the price around based on global market expectations.  The second bath is the fractional Brownian motion which is essentially the individual users interacting with the price in a way that is inspired by their memory.  Hence, our full fractional stochastic differential equation is

$$
\lambda \frac{d^2 x}{dt^2} + \beta \lambda \frac{dx}{dt} - a D^{\alpha} \big[ x(t) \big] = \eta(t) - \xi_{\alpha} (t),
\label{eq:fin-lang-eq}
$$

An interesting side-note would be that the values for $\Theta_1$ and $\Theta_2$ are given by the fluctuation-dissipation theorem and are based on the financial system responding to noise in an expected manner.  It would be interesting to investigate an out-of-equilibrium system when those values for $\Theta_1$ and $\Theta_2$ do not hold, especially given the work in Robin Verstraten's master's thesis on out-of-equilibrium systems.

### Financial Intuition from the Limiting case of No Memory
Now, instead of taking the kernel $K_R$ to obtain the Caputo fractional derivative, we consider that the kernel has no memory, then the fractional derivative goes to integer order and the differential equation becomes

$$
\lambda \frac{d^2 x}{dt^2} + (\beta \lambda - a) \frac{dx}{dt} = f(t),
$$

where according to the fluctuation-dissipation theorem and stable solution existence, $\beta \lambda \geq a$.  The left-hand side is the rate of change of demand when the price moves a single unit which is a dynamical variable describing the entire market.  This describes how the demand is expected to change in time when the price moves a single unit based purely on the market depth.  When the relationship is strictly equal, it would mean that the impact of recent trends on user's actions exactly matches up with how the market expects users to respond.  This is an entirely capitalistic model in which users exclusively seek profits.  The friction term goes to zero and we are left with a simple stochastic differential equation, but this is not interesting for us.

When $\beta \lambda < a$, the impact of recent trends is less than what the market expects.  An idealistic trader in the market will hold onto the asset independent of how the price is behaving because they are a core believer in whatever the asset represents.  A timely example of this would be Bitcoin: there exists Bitcoin users who believe so strongly that it is the future so no matter what the price is, they will never sell.  The impact of recent trends is essentially averaged over all users in the market so if there are idealists in the market, they do no contribute to the parameter $a$ and thus the value is less than what the market depth would predict.  Therefore, the existence of idealists in the market is what generates friction in this model.  The term used here is idealists, but that implies an assumption of their motivation.  These actors can also been seen as non-optimal users in the market, i.e., a kind of disorder.

## Monte Carlo Solution Finder
At this point, we have a differential equation that describes a financial system with memory and a numerical method for generating fractional Brownian motion.  Now, we solve the equation.  The financial differential equation is quite difficult to solve analytically, so instead we decide to solve it numerically.  The key to most numerical equation solvers is how the next guess at the solution is made.  The method of choosing the next guess is typically built around what kind of equation is trying to be solved.  However, this can require extensive knowledge of the equation being solved.  For some ease in building a numerical solver, we use a simple Metropolis-Hastings Monte Carlo method to find our next guess.  The method goes as follows,

include figure of mc pseudo-code

The function $g(t)$ is the full differential equation without the fractional derivative term.  Thus, it is straight-forward to find an analytic solution.  The interesting part is to figure out how to guess the next possible solution.  This is where the Monte Carlo method comes in.

The residue for each time point is taken to effectively be the energy of that point.  A Monte Carlo method works by minimizing the energy.  The algorithm proposes a random change in the price position at that time point, with a given hard-coded maximum change, and calculates the new residues for all points in time.  If the new price position generates a residue that is smaller than those of the previous time point, then the new price position is accepted.

Normally, a Monte Carlo code would attempt to move every single time point at each Monte Carlo step, however our code does not do this.  The residuals calculate how close a given time point is to the true solution and because they depend on history, the residuals of a given time point depend on the time points in the past.  If the time points in the past are not within the prescribed tolerance to the true solution, then the residuals calculated do not accurately represent the closeness of this time point to its true solution.  An attempt to move this time point would then be useless because its acceptance does not guarantee that it is now closer to the solution.

To avoid this issue, moves are only attempted for a small group of time points, as seen in the \textbf{Next Guess Pseudo-code}.  In this way, the solution is found by effectively moving through time so that the history used to calculate the residuals is actually representative of the true history of a given time point.

## Results
In order to analyze the behavior of our financial model we examine the Mean Squared Displacement (MSD),

$$
\textrm{MSD}(t) = \big\langle \vert x(t) - x(t=0) \vert^2 \big\rangle,
$$

where $\langle...\rangle$ is an ensemble average over a given number of different paths.  The MSD is effectively equivalent to the variance and thus can be used to calculate volatility and to examine the dispersion relations of our system.

Unless otherwise specified, the system begins with an initial velocity which is related to the thermal energy in the system.  This is done in order to start-off the system in thermal equilibrium because we are assuming that the system always follows the equipartition theorem.  Therefore, the initial velocity is always $v_0=\sqrt{k_B T/\lambda}$.

### Liquidity and Volatility
As previously discussed, the liquidity of a given financial market plays an important role in the volatility of that market.  Therefore it is important that our model reproduces these same relationship.  One definition of volatility is the variance over a finite period of time, and is because the MSD is equivalent to the variance, the MSD can be used to calculate volatility.

volat-liquid figure

Figure volat-liquid shows the 100-time-step volatility as a function the market liquidity parameter, $\lambda$, for a few values of the fractional derivative order $\alpha$ where all other parameters are kept the same.  As expected, as the market liquidity increases, the finite time volatility decreases showing that more liquid markets are more stable than less liquid markets.  The error bars here are quite large because the volatility comes from the MSD which is only averaged over 20 paths.  This comes from the computational limits of our simulation.  However, the trend between liquidity and volatility is still clear enough to say the relationship does exist.

%There are two types of volatile markets.  The first type is a market in which the price moves in large jumps that always go the same direction and thus has high volatility.  The magnitude of the MSD in a market of this type will be large.  The inverse is also true for low volatility markets of this kind.  The second type is where large jumps take place but they happen in both directions and so the MSD is quite small because these jumps cancel each other out.  Both of these market types are described by our model and so the calculated volatility will represent both types.  The figure shows that as the liquidity, $\lambda$, increases, the magnitude of the MSD, which is the 100-time-step volatility, decreases as expected.  When liquidity is low, this means that type one markets will have high volatility, and type two markets will have low volatility.  Therefore, because our model takes into account both types, the error bars of the volatility are quite large because they are incorporating both market types which have very different volatilities.  As the liquidity increases, the type one market's volatility will decrease, getting closer to the volatility of type two markets and so the error bars get smaller.

However, the large error bars suggest that is no relationship between the volatility and the fractional derivative order.  Indeed, the existence of memory in the system does not play a significant enough role to change the fundamental relationship between volatility and liquidity.

### Expected MSD Behavior
For all MSD plots, the parameters being used are $\lambda=500$, $\beta=1.0=a$ and the time scale is in seconds. To understand the asymptotic behavior of the financial differential equation, we look at four different interaction components: fractional derivative with colored noise, fractional derivative with white noise, integer derivative with colored noise, and integer derivative with white noise.  We know that the fractional derivative with colored noise and the integer derivative with white noise are governed by the fluctuation-dissipation theorem.

In the limit of no memory, the fractional derivative with colored noise can be neglected, and the regular integer order Langevin equation is retrieved, which leads to regular Brownian motion.  The noise term on the right-hand side is white noise and the friction term is a first order derivative.

From Lutz, when the Langevin equation contains just a fractional derivative and colored noise term, the scaling goes as $t^{\alpha}$ in the long-term.  The limiting case for this where $\alpha \rightarrow 1$ also maintains the linear scaling of regular Brownian motion.

However, because our system contains both a fractional derivative and an integer derivative, there exists more complex interactions between these terms. From the master's thesis of Robin Verstraten we also know what happens with a fractional derivative interacting with a white noise term.  The interaction is out-of-equilibrium behavior which does not follow the fluctuation-dissipation theorem, but instead has meta-stable states as a function of the fractional derivative order, $\alpha$.  The long-term scaling of the MSD then goes as $t^{2\alpha-1}$ which is why it is shown for comparison in the table of all MSDs.

### Dispersion Regions

figure msd 0p2

The figure above shows the MSD for a fractional derivative order of $\alpha=0.2$, which has a positive correlation to history, and a finite initial velocity.  The figure shows three distinct regions each with different dispersion relations.  At the start, the curve goes as $t^2$ which is ballistic dispersion.  These are the times where the trace "particle" is exploring the environment and does not yet have any memory.

The third region, at the end, is the long-term behavior of the system.  In the case of $\alpha=0.2$, the dispersion goes as $t^{2-\alpha}$.  A discussion of the properties of this phase will follow.  Between the beginning exploratory behavior and the long-term behavior, there is a transition region.  Here, the MSD plateaus becoming a constant which implies more solid-like behavior.  We believe this is where the system is building up more and more memory.  The memory is still quite small, thus it has very little impact on the dynamics.  Once the memory has built up to a significant enough level, the system moves into the third region which displays the dynamics of a system with memory.  The MSD figure with no memorry shows that for a system with no memory and there is no middle transition region.  Therefore, the existence of a transition region appears tied to the inclusion of memory in the system.  For regular Brownian motion, the scaling relation with time is just linear and this is exactly the long-term behavior seen in the figure.

figure no memory

From a financial perspective, this could mean that there is a reaction time and shape associated to a given financial system.  

### Long-term Phase Behavior

figure 0.3 and 0.4

The figure above shows the MSD curves for two values of $\alpha$ both with finite initial velocity.  In the third region, which corresponds to the long-term phase behavior, the two curves follow different dispersion relations.  For $\alpha=0.3$, it scales as $t^{2-\alpha}$ while for $\alpha=0.4$ it scales as $t^{\alpha}$.  It appears that around these values of $\alpha$ the dynamics of the dispersion change into a different form of anomalous diffusion.  In Figure \ref{fig:msd-phases-trans}, $\alpha=0.4$ does look close to linear scaling, but the MSD table shows the full range of $\alpha$ parameters.  After $\alpha=0.4$, it can be seen that the MSD scales as $t^{\alpha}$ for the rest of the $\alpha$ values.

#### Small $\alpha$ Behavior
For small values of $\alpha$, as in the $\alpha=0.2$ figure, there appears to to be kind of staircase shape to the long-term MSD.  Since the volatility can be calculated from the MSD, it appears that the volatility goes through repeating regions of flat volatility, then increasing volatility.  While the scaling of the long-term behavior shows anomalous diffusion, the overall staircase shape is similar to that of a marginal glass.  This region is not exactly a marginal glass, but some of its qualitative features are similar and so can be used to potentially infer what is physically taking place in the system.

The MSD appears to move from flatter regions, which are more glass-like behavior, to higher volatility regions which correspond to anomalous diffusion, and then back to glass-like again.  In a marginal glass, this behavior takes place because single particles cannot move on their own, but they become part of a group of particles that can move as a whole.  These cages have a fractal structure with cages inside of cages.  Because of this fractal structure, there exist metastable states for the particles given by the size of the cages.  When a region of anomalous diffusion reaches the next cage size, the MSD saturates and plateaus appear.  However, because the states are metastable, the plateau only exists for a finite period of time.  More study is required to see which parameters impact the lifetime of the plateaus and the cage sizes.  The existence of this qualitative behavior in the MSD shows that the topic is worth further investigation.  In order to properly characterize this phase, the asymptotic limits need to be understood, which would require more computational power.  A marginal glass is not stable in the long term because that would require infinite energy so these asymptotic limits are important.

It is difficult to say what these cages represent in a financial system as the system is very complex and we will not take a random guess with sufficient reasoning behind the analysis.  However, if these regions can be understood and recognized, it could lead to a better understanding of changes in market dynamics.

figure 0.18 and 0.6


\subsection{Oscillatory Behavior}
As the fractional derivative order increases, the strength of the memory increases as well; memory becomes more important.  The figure above shows MSD plots for $\alpha=0.18$ and $0.6$ with all other parameters kept the same.  For larger values of $\alpha$, there appears to be larger oscillations in the MSD.  This leads us to believe that, in a rudimentary way, it is the memory that induces oscillations in the MSD.

There is not enough precision is our calculations to say whether these oscillations are truly periodic or if they are random.  However, the existence of large changes in the MSD can be seen as changes in the volatility over time.  Volatility is in part a reflection of the confidence users have in a given financial market; high volatility means users are unsure about the future while low volatility means high confidence in future expectations.  Because the volatility calculated in this way comes from historical information about the asset price, it is the realized volatility.

Another form of volatility is the implied volatility which can be calculated from options pricing models.  One of the most famous options pricing models is the Black-Scholes model which is based on a stochastic process that represents regular Brownian motion.  Since the development of fractional Brownian motion, many options pricing models have been built using this method and are typically called rough volatility models.  These include the rough Heston model, the SABR model, the rough Bergomi model, along with many others.

The mathematical base of these models are stochastic differential equations, however the model described in this report has a different structure.  It is possible that our model is equivalent to one of the other stochastic models already developed, but it is a non-trivial task to show that equivalence.  By this same reasoning, it is also non-trivial to use the mathematical tools of other fractional volatility models to analyze the results found from our model.

### Asymptotic Behavior
Given our expectations for the scaling of the MSD, we now analyze the simulation results more closely.  For $\alpha = 0.6$, the scaling appears to be closer to $t^{\alpha}$, however, as $\alpha$ grows the scaling starts to look more like $t^{2\alpha-1}$, for example $\alpha = 0.9$.  When $\alpha > 0.5$, it appears that the fractional derivative term is dominant where low values correspond to more interaction with the colored noise and higher values connect to the white noise.

For $\alpha = 0.2$ and $0.3$, the scaling is clearly $t^{2-\alpha}$, while as $\alpha$ increases the scaling gets closer to $t^{\alpha}$.  This leads us to believe that, qualitatively, the integer derivative interaction with colored noise leads to the $t^{2-\alpha}$ behavior.

The dispersion relation does transition between these different phases.  However the path by which it does this is unclear.  The complexity of having many terms in the differential equation means that it is difficult to analytically describe the transition behavior.  More extensive research would be required to analyze the equation and understand the expected asymptotic behavior and its transition as a function of $\alpha$.

### Conclusions
The model reproduces the fundamental relationship between volatility and liquidity, as seen in the figure.  This relationship is unaffected by the existence of memory within the system.  For small values of the fractional derivative order, there exists behavior qualitatively similar to a marginal glass.  While this is far from conclusive, it shows some interesting features of financial markets that warrant further study.  Larger values of the fractional derivative order lead to oscillations in the realized volatility.  Because of the mathematical structure of the model it is difficult to use the typical tools of rough volatility financial models.  Therefore, more research would be required to compare this model with those already in existence to see if there are truly any different results.  The model analyzed in this report has generated potentially interesting behavior, but it has raised more questions than answers: continued study into this model is necessary to understand better what it is describing.
