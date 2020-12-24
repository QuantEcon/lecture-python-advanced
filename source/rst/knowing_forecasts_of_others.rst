.. _knowing_the_forecast_of_others_v3:

.. include:: /_static/includes/header.raw

.. highlight:: python3


**************************************
Knowing the Forecasts of Others
**************************************


.. contents:: :depth: 2


:cite:`lucas75`, :cite:`kasa`, and
:cite:`townsend` demonstrated that the assumption that
decision makers have incentives to infer hidden persistent state
variables from equilibrium prices and quantities is a potential source
both for additional impulses and for elongated impulse response
functions in business cycle models. :cite:`townsend`
indicated that models that incorporate such incentives can naturally
induce decision makers in effect to forecast the forecast of others.
This theme has been pursued and extended in recent analyses in which
decision maker’s imperfect information forces them into pursuing an
infinite recursion of forming beliefs about the beliefs of other
(e.g., :cite:`ams`).

:cite:`lucas75` side stepped the problem of forecasting the
forecasts of others by assuming that decision makers pool their
information before forecasting. Because he didn’t assume such pooling,
confronted the forecasting the forecasts of others problem. However,
that led to what he thought was an intractable state space. That led him
to proposed more manageable model that he argued could do a good job of
approximating the intractable model.

By applying technical machinery of :cite:`PCL`,
:cite:`PS2005` showed that there is a recursive
representation of the equilibrium of the perpetually and symmetrically
uninformed model formulated but not completely solved in section 8 of
:cite:`townsend`. Their computational method is recursive:
it combines the Kalman filter with invariant subspace methods for
solving systems of Euler
equations [#footnote1]_ . As :cite:`singleton`,
:cite:`kasa`, and :cite:`sargent91` also
found, the equilibrium is fully revealing: observed prices tell
participants in industry :math:`i` all of the information held by
participants in market :math:`-i` (:math:`-i` means not :math:`i`). This
means that higher-order beliefs play no role: seeing equilibrium prices
in effect lets decision makers pool their information
sets [#footnote2]_ . The disappearance of higher order beliefs means that
decision makers in this model do not really face a problem of
forecasting the forecasts of others. They know those forecasts because
they are the same as their own.

The presence of a common hidden state variable is the only thing that
inspires decision makers in one market to condition their decisions on
the history of prices in the other market.

:cite:`townsend` noted that in his model with perpetually
and symmetrically uninformed decision makers, the dimension of the state
space seemed to explode because it seemed to be necessary for decision
makers to keep track of an infinite history of vectors of observables.

That *curse of dimensionality* deterred Townsend from characterizing or
computing an equilibrium of that model.

Instead he constructed another model and computed its equilibrium. To
construct this model he assumed that after a finite number :math:`j`
periods, the (lagged) value of the key hidden state variable is revealed
to the decision maker.

:cite:`sargent91` proposed a way to compute an equilibrium
without making Townsend’s approximation. Extending the reasoning of
:cite:`muth1960`, Sargent noticed that it is possible to
summarize the relevant history with a low dimensional object, namely, a
small number of current and lagged forecasting errors. Positing an
equilibrium in a space of perceived laws of motion for endogenous
variables that takes the form of a vector autoregressive, moving
average, Sargent described an equilibrium as a fixed point of a mapping
from the perceived law of motion to the actual law of motion of that
form. Sargent worked in the time domain and had to guess and verify the
appropriate orders of the autoregressive and moving average pieces of
the equilibrium representation. However, by working in the frequency
domain :cite:`kasa` showed how to discover the appropriate
orders of the autoregressive and moving average parts, and also how to
compute an equilibrium.

Our recursive computational method, which stays in the time domain, also
discovers the appropriate orders of the autoregressive and moving
average pieces. In addition, by displaying equilibrium representations
in the form of :cite:`PCL`, :cite:`PS2005`
showed how the moving average piece is linked to the innovation process
of the hidden persistent component of the demand shock. That scalar
innovation process is the additional state variable contributed by the
problem of extracting a signal from equilibrium prices that decision
makers face in Townsend’s model.

**Tom**: CONSOLIDATE AND EDIT ABOVE AND BRING SOME OF IT TO CONCLUDING
SECTION

This lecture describes a model of two industries that are linked in a
single way: shocks to the demand curves for their products have a common
component.

We shall decompose the model into components that we shall then put
together in several ways that help us appreciate the structure of the
**pooling equilibrium** model that ultimately concerns us.

While keeping other aspects of the model the same, we shall to study
consequences of alternative assumptions about what decision makers
observe.

The Setting
============

Firms in each of two industries :math:`i=1,2` use a single factor of
production, capital :math:`k_t^i`, to produce output of a single good,
:math:`y_t^i`.

Firms bear quadratic costs of adjusting their capital stocks.

We let capital letters denote market wide objects and lower case letters
denote objects chosen by a representative firm.

To rationalize the big :math:`K`, little :math:`k` connection, we can
think of there being a continua of each type of firm, each indexed by
:math:`\omega \in [0,1]` with
:math:`K^i = \int_0^1 k^i(\omega) d \omega`.

A representative firm in industry :math:`i` has production function
:math:`y_t^i = f k_t^i`, :math:`f >0`, acts as a price taker with
respect to output price :math:`P_t^i`, and maximizes

.. math::
    :label: town1

      \begin{aligned}
      E_0^i \sum_{t=0}^\infty \beta^t \left\{ P_t^i f k_t^i - .5
         h (k_{t+1}^i - k_t^i)^2 \right\} ,
      \quad h >0  .\end{aligned}

Demand in industry :math:`i` is described by the inverse demand curve

.. math::
    :label: town2

      \begin{aligned}
      P_t^i = - b Y_t^i + \theta_t + \epsilon_t^i  , \quad b >0,
      \end{aligned}

where :math:`P_t^i` is the price of good :math:`i` at :math:`t`,
:math:`Y_t^i = f K_t^i` is output in market :math:`i`, :math:`\theta_t`
is a persistent component of a demand shock that is common across the
two industries, and :math:`\epsilon_t^i` is an industry specific
component of the demand shock that is i.i.d. and whose time :math:`t`
marginal distributon is :math:`{\mathcal N}(0, \sigma_{\epsilon}^2)`.

We assume that :math:`\theta_t` is governed by

.. math::
    :label: town2a

      \begin{aligned}
      \theta_{t+1} = \rho \theta_t + v_{t}
      \end{aligned}

where :math:`\{v_{t}\}` is an i.i.d. sequence of Gaussian shocks each
with mean zero and variance :math:`\sigma_v^2`.

To simplify notation, we’ll study a special case of the model by setting
:math:`h=f=1`.

Big K, little k connection
-----------------------------

In equilibrium, :math:`k_t^i = K_t^i`, but as usual we must distinguish
between :math:`k_t^i` and :math:`K_t^i` when we pose the firm’s
optimization problem.

Following Townsend, we eventually want to assume that at time :math:`t`
firms in industry :math:`i` observe
:math:`k_t^i, Y_t^i, P_t^i, (P^{-i})^t`, where :math:`(P^{-i})^t` is the
history of prices in the other market up to time :math:`t`.

Because the representative firm :math:`i` sees the price as well as the
aggregate state variable :math:`Y_t^i` in its own industry, it can infer
the total demand shock :math:`\theta_t + \epsilon_{t}^i`.

However, at time :math:`t`, the firm sees only :math:`P_t^{-i}` and does
not see :math:`Y_t^{-i}`, so that firm :math:`i` does not directly
observe :math:`\theta_t + \epsilon_t^{-i}`.

Punch line
============

Nevertheless it turns out that in the **pooling equilibrium** that
ultimately interests us, a firm in industry :math:`i` will be able to
infer the composite shock :math:`\theta_t + \epsilon_t^{-i}` from the
history of random variables that it observes at :math:`t`.

We shall proceed to establish this result and others in steps.

Strategy
-----------

To prepare to compute a pooling equilibrium, we shall first compute
equilibrium laws of motion for capital in industry :math:`i` under a
sequence of assumptions about what a representative firm observes.

Successive members of this sequence make a representative firm’s
information more and more obscure.

We begin with the most information, then gradually withdraw information
in a way that approaches and eventually reaches the information
structure that that we are ultimately interested in.

Thus, we shall compute equilibria under the following alternative
information structures:

-  **Perfect foresight:** future values of
   :math:`\theta_t, \epsilon_{t}^i` are observed in industry :math:`i`.

-  **Observed but stochastic** :math:`\theta_t`:
   :math:`\{\theta_t, \epsilon_{t}^i\}` are realizations from a
   stochastic process; current and past values of each are observed at
   time :math:`t` but future values are not.

-  **One noise-ridden observation on** :math:`\theta_t`: Values of
   :math:`\{\theta_t, \epsilon_{t}^i\}` are never observed. However, at
   time :math:`t`, a history :math:`w^t` of a scalar noise-ridden
   observations on :math:`\theta_t` is observed at time :math:`t`.

-  **Two noise-ridden observations on** :math:`\theta_t`: Values of
   :math:`\{\theta_t, \epsilon_{t}^i\}` are never observed. However, at
   time :math:`t`, a history :math:`w^t` of *two* noise-ridden
   observations on :math:`\theta_t` is observed at time :math:`t`.

Successive computations build one on another.

We proceed by first finding an equilibrium under perfect foresight.

To compute an equilibrium with :math:`\theta_t` observed, we use a
*certainty equivalence principle* to justify modifying the perfect
foresight equilibrium by replacing future values of
:math:`\theta_s, \epsilon_{s}^i, s \geq t` with mathematical
expectations conditioned on :math:`\theta_t`.

This provides the equilibrium when :math:`\theta_t` is observed at
:math:`t` but future :math:`\theta_{t+j}`\ s and
:math:`\epsilon_{t+j}^i`\ s are not observed.

To find an equilibrium when only a history :math:`w_t` of a single noise
ridden observations on :math:`\theta_t` is observed, we again apply a
certainty equivalence principle and replace future values of the random
variables :math:`\theta_s, \epsilon_{s}^i, s \geq t` with their
mathematical expectations conditioned on :math:`w^t`.

To find an equilibrium when only a history :math:`w_t` of a *two* noisy
signal on :math:`\theta_t` is observed, we replace future values of the
random variables :math:`\theta_s, \epsilon_{s}^i, s \geq t` with their
mathematical expectations conditioned on history :math:`w^t`.

The equilibrium with two noise-ridden observations on :math:`\theta_t`
we call a **pooling equilibrium**

-  It corresponds to an arrangement in which at the beginning of each
   period firms in industries :math:`1` and :math:`2` somehow get
   together and share information about current values of their noisy
   signals on :math:`\theta`.

We want ultimately to compare outcomes in such a *pooling equilibrium*
with an equilibrium under the following information structure for a firm
in industry :math:`i` that interested Townsend:

-  **Firm** :math:`i`\ ’s **noise-ridden signal on** :math:`\theta_t` **and the
   price in industry** :math:`-i`**:** At time :math:`t`, a firm in industry
   :math:`i` observes a history :math:`w^t` of *one* noise-ridden signal
   on :math:`\theta_t` and a history of industry :math:`-i`\ ’s price is
   observed.

XXXX **Tom revise and move stuff below** match the equilibrium that we
seek because equilibrium prices in that equilibrium completely reveal to
firms in industry :math:`i` the noisy signal about the demand shock
received by firms in industry :math:`-i`.

In this way, we construct benchmarks against which we can interpret an
equilibrium to under the following informatin structure:

Townsend’s model XXX that we shall compute in section
:eq:`PCL2` by applying the machinery of
:cite:`PCL`.

Equilibrium conditions
=======================

It is convenient to formulate the firm’s problem as a discrete time
Hamiltonian by forming the Lagrangian for the problem without
uncertainty:

.. math::

   \begin{aligned}
   J=\sum_{t=0}^\infty \beta^t \left\{
   P_t^i  k_t^i - .5   (\mu_t^i)^2  + \phi_t^i \left[
      k_t^i + \mu_t^i - k_{t+1}^i \right]  \right\} \end{aligned}

where :math:`\{\phi_t^i\}` is a sequence of Lagrange multipliers on the
transition law for :math:`k_{t+1}^i`. First order conditions for the
nonstochastic problem are

.. math::
    :label: town4

      \begin{aligned}
      \phi_t^i & =  \beta \phi_{t+1}^i + \beta  P_{t+1}^i   \\
      \mu_t^i & =  \phi_t^i .   \end{aligned}


Substituting the demand function :eq:`town2` for
:math:`P_t^i`, imposing the condition that the representative firm is
representative ( :math:`k_t^i = K_t^i`), and using the definition below
of :math:`g_t^i`, the Euler equation :eq:`town4`, lagged
by one period, can be expressed as
:math:`- b k_t^i + \theta_t + \epsilon_t^i + (k_{t+1}^i - k_t^i) - g_t^i =0`
or

.. math::
    :label: pcl11


      \begin{aligned}
      k_{t+1}^i = (b+1) k_t^i - \theta_t - \epsilon_t^i + g_t^i
      \end{aligned}

where we define :math:`g_t^i` by

.. math::
    :label: town7

      g_t^i = \beta^{-1} (k_t^i - k_{t-1}^i)
      
We can write the Euler equation :eq:`town4` in terms of :math:`g_t^i`:


.. math::
    :label: pcl10

      \begin{aligned}
      g_t^i = P_t^i + \beta g_{t+1}^i
      \end{aligned}

In addition, we have the law of motion for :math:`\theta_t`,
:eq:`town2a`, and the demand equation
:eq:`town2`.

In summary, with perfect foresight, equilibrium conditions for industry
:math:`i` consist of the following system of difference equations:


    

.. math::
    :label: sol1

      \begin{aligned}
      k_{t+1}^i & =  (1+b)k_t^i - \epsilon_t^i -\theta_t + g_t^i \\
      \theta_{t+1} & =  \rho \theta_t + v_t \\
      g_{t+1}^i  & = \beta^{-1} (g_t^i - P_t^i)   \\
      P_t^i & =  -b k_t^i + \epsilon_t^i + \theta_t  \end{aligned}

Without perfect foresight, the same system prevails except that the
following equation replaces the third equation of :eq:`sol1`:

.. math::

   \begin{aligned}
   g_{t+1,t}^i = \beta^{-1} (g_t^i - P_t^i) \end{aligned}


where
:math:`x_{t+1,t}` denotes the mathematical expectation of
:math:`x_{t+1}` conditional on information at time :math:`t`.

Solution under perfect foresight
------------------------------------

Our first step is to compute the equilibrium law of motion for
:math:`k_t^i` under perfect foresight.

Let :math:`L` be the lag
operator. [#footnote3]_

Equations :eq:`pcl10` and :eq:`pcl11`
imply the second order difference equation in
:math:`k_t^i` [#footnote4]_ .


.. math::
    :label: euler1

      \begin{aligned}
      \left[ (L^{-1} - (1+b))(1-\beta L^{-1}) + b\right] k_t^i
      = \beta L^{-1} \epsilon_t^i + \beta L^{-1} \theta_t .
      \end{aligned}

Factor the polynomial in :math:`L` on the left side as:

.. math::

      \begin{aligned}
      -\beta [L^{-2} -(\beta^{-1} + (1+b))L^{-1} + \beta^{-1}]
      = \tilde \lambda^{-1}(L^{-1} - \tilde \lambda)(1-\tilde \lambda \beta L^{-1})
      \end{aligned}

where :math:`|\tilde \lambda | < 1` is the smaller root and
:math:`\lambda` is the larger root of
:math:`(\lambda-1)(\lambda-1/\beta)=b\lambda`.

Therefore, :eq:`euler1` can be expressed as

.. math::

      \begin{aligned}
      \tilde \lambda^{-1}(L^{-1} - \tilde \lambda) (1-\tilde \lambda \beta L^{-1})
      k_t^i = \beta L^{-1} \epsilon_t^i + \beta L^{-1} \theta_t .
      \end{aligned}

Solving the stable root backwards and the unstable root forwards gives

.. math::

      \begin{aligned}
      k_{t+1}^i = \tilde \lambda k_t^i + {\tilde \lambda \beta \over 1 -\tilde
      \lambda \beta L^{-1}}
      (\epsilon_{t+1}^i + \theta_{t+1}  )
      \end{aligned}

Thus, under perfect foresight the capital stock satisfies

.. math::
    :label: town5

      \begin{aligned}
      k_{t+1}^i = \tilde \lambda k_t^i + \sum_{j=1}^\infty (\tilde \lambda \beta)^j
      (\epsilon_{t+j}^i +  \theta_{t+j}) .
      \end{aligned}

Next, we shall use alternative forecasting formulae in
:eq:`town5` to compute the equilibrium decision rule
under alternative assumptions about the information available to
decision makers in market :math:`i`.

Solution with :math:`\theta_t` stochastic but observed at :math:`t`
======================================================================

If future :math:`\theta`\ ’s are unknown at :math:`t`, it is appropriate
to replace all random variables on the right side of
:eq:`town5` with their conditional expectations based on
the information available to decision makers in market :math:`i`.

For now, we assume that this information set
:math:`I_t^p = \begin{bmatrix} \theta^t & \epsilon^{it} \end{bmatrix}`,
where :math:`z^t` represents the infinite history of variable
:math:`z_s` up to time :math:`t`.

Later we shall give firms less information about :math:`\theta_t`.

To obtain the counterpart to :eq:`town5` under our
current assumption about information, we apply a certainty equivalence
principle.

In particular, it is legitimate to take :eq:`town5` and
replace each term :math:`( \epsilon_{t+j}^i+ \theta_{t+j} )` on the
right side with
:math:`E[ (\epsilon_{t+j}^i+ \theta_{t+j}) \vert \theta^t ]`.

After using :eq:`town2a` and the i.i.d. assumption about
:math:`\{\epsilon_t^i\}`, this gives

.. math::

   \begin{aligned}
   k_{t+1}^i = \tilde \lambda k_t^i + {\tilde \lambda \beta \rho \over 1 -
   \tilde \lambda \beta \rho }
   \theta_t
   \end{aligned}

or

.. math::
    :label: solution1

      \begin{aligned}
      k_{t+1}^i = \tilde \lambda k_t^i  + {\rho \over  \lambda - \rho} \theta_t
      \end{aligned}

where :math:`\lambda \equiv (\beta \tilde \lambda)^{-1}`.

For future purposes, it is useful to represent the solution for
:math:`k_t^i` recursively as


.. math::
    :label: sol10

      \begin{aligned}
      k_{t+1}^i  & =  \tilde \lambda k_t^i  + {1 \over \lambda - \rho}
      \hat \theta_{t+1}  \\
      \hat \theta_{t+1}  & =  \rho \theta_t \\
      \theta_{t+1} & =  \rho \theta_t + v_t.   \end{aligned}

Filtering
-----------

One noisy signal
^^^^^^^^^^^^^^^^^^^^^^^^^^

We get closer to the model that we ultimately want to study by now
assuming that firms in market :math:`i` do not observe :math:`\theta_t`,
but instead observe a history of noisy signals :math:`w^t`.

In particular, assume that

.. math::
    :label: kf1&2

      \begin{aligned}
      w_t  & =   \theta_t + e_t  \label{kf1}  \\
      \theta_{t+1} & =  \rho \theta_t + v_t 
      \end{aligned}

where :math:`e_t` and :math:`v_t` are mutually independent
i.i.d. Gaussian shock processes with means of zero and variances
:math:`\sigma_e^2` and :math:`\sigma_v^2`, respectively.

Define

.. math::

   \begin{aligned}
   \hat \theta_{t+1} = E(\theta_{t+1} | w^t)
   \end{aligned}

where :math:`w^t` denotes the history of the :math:`w_s` process up to
and including :math:`t`.

Associated with the state-space representation
:eq:`kf1&2` is the *innovations
representation*


.. math::
    :label: kf3&4

      \begin{aligned}
      \hat \theta_{t+1}  & =     \rho \hat \theta_t + k a_t  \\
      w_t & =  \hat \theta_t + a_t 
      \end{aligned}

where :math:`a_t \equiv w_t - E(w_t | w^{t-1})` is the *innovations*
process in :math:`w_t` and the Kalman gain :math:`k` is

.. math::
    :label: kal1

      \begin{aligned}
      k = {\rho p \over p + \sigma_e^2} \end{aligned}

and where :math:`p` satisfies the Riccati equation

.. math::
    :label: kf6

      \begin{aligned}
      p = \sigma_v^2   + { p \rho^2 \sigma_e^2 \over \sigma_e^2 + p}.
      \end{aligned}

Define the state *reconstruction error* :math:`\tilde \theta_t` by

.. math::
   
   \begin{aligned}
   \tilde \theta_t = \theta_t - \hat \theta_t .
   \end{aligned}

Then :math:`p = E \tilde \theta_t^2`. Equations :eq:`kf1&2`
and :eq:`kf3&4` imply

.. math::
    :label: kf7

      \begin{aligned}
      \tilde \theta_{t+1} = (\rho - k) \tilde \theta_t + v_t - k e_t .
      \end{aligned}

Now notice that we can express :math:`\hat \theta_{t+1}` as


.. math::
    :label: kf8

     \hat \theta_{t+1}   = [\rho \theta_t + v_t]  + [ ke_t - (\rho -k) \tilde \theta_t - v_t]  ,
      

where the first term in braces  equals
:math:`\theta_{t+1}` and the second term in braces equals
:math:`-\tilde \theta_{t+1}`.

Additional state variable: :math:`\theta`-reconstruction error:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can express :eq:`solution1` as


.. math::
    :label: solution2

      \begin{aligned}
      k_{t+1}^i = \tilde \lambda k_t^i + {1 \over \lambda - \rho}
      E \theta_{t+1} | \theta^t .
      \end{aligned}

An application of a certainty equivalence principle asserts that when
only :math:`w^t` is observed, the appropriate solution is found by
replacing the information set :math:`\theta^t` with :math:`w^t` in
:eq:`solution2`.

Making this substitution and using :eq:`kf8` leads to

.. math::
    :label: kf9

      \begin{aligned}   
      k_{t+1}^i   = \tilde \lambda k_t^i + {\rho \over  \lambda - \rho} \theta_t + {k \over  \lambda - \rho} e_t  - {\rho - k \over  \lambda - \rho} \tilde \theta_t .
      \end{aligned}
      

Simplifying equation :eq:`kf8`, we also have

.. math::
    :label: kf8a

      \begin{aligned}
      \hat \theta_{t+1}  = \rho \theta_t +  ke_t - (\rho -k) \tilde \theta_t  .
      \end{aligned}
      

Equations :eq:`kf9`, :eq:`kf8a` describe
the solution when :math:`w^t` is observed.

Relative to :eq:`solution1`, the solution acquires a new
state variable, namely, the :math:`\theta`–reconstruction error,
:math:`\tilde \theta_t`.

For future purposes, by using :eq:`kal1`, it is useful to
write :eq:`kf9` as

.. math::
    :label: sol2a


      \begin{aligned}
      k_{t+1}^i = \tilde \lambda k_t^i + {\rho \over  \lambda - \rho } \theta_t + {1 \over  \lambda - \rho} {p \rho \over p + \sigma_e^2} e_t - {1 \over \lambda - \rho} {\rho \sigma_e^2 \over p + \sigma_e^2}
      \tilde \theta_t
      \end{aligned}

In summary, when decision makers in market :math:`i` observe a noisy
signal :math:`w_t` on :math:`\theta_t` at :math:`t`, we can represent an
equilibrium law of motion for :math:`k_t^i` as

.. math::

      \begin{aligned}
      k_{t+1}^i & =  \tilde \lambda k_t^i + {1 \over \lambda - \rho}
      \hat \theta_{t+1} \label{sol4;a} \\
      \hat \theta_{t+1} & =  \rho \theta_t + {\rho p \over p + \sigma_e^2} e_t
         - {\rho \sigma_e^2 \over p + \sigma_e^2} \tilde \theta_t
      \label{sol4;b} \\
      \tilde \theta_{t+1} & = { \rho \sigma_e^2 \over p + \sigma_e^2} \tilde
         \theta_t - {p \rho \over p + \sigma_e^2} e_t + v_t
      \label{sol4;c} \\
      \theta_{t+1} & =  \rho \theta_t + v_t .  \label{sol4;d} \end{aligned}

Two noisy signals
--------------------

We now construct a **pooling equilibrium** by assuming that a firm in
industry :math:`i` receives a vector :math:`w_t` of *two* noisy signals
on :math:`\theta_t`:


.. math::

      \begin{aligned}
      \theta_{t+1} & =  \rho \theta_t + v_t  \label{kf20} \\
      w_t   & = 
      \begin{bmatrix} 1 \\ 1 \end{bmatrix}
      \theta_t
         + \begin{bmatrix} e_{1t} \\ e_{2t} \label{kf21} \end{bmatrix}
      \end{aligned}

To justify that we are constructing is a **pooling equilibrium** we can
assume that

.. math::

   \begin{aligned}
   \begin{bmatrix} e_{1t} \\ e_{2t}  \end{bmatrix} =
   \begin{bmatrix} \epsilon_{t}^1 \\ \epsilon_{t}^2  \end{bmatrix}
   \end{aligned}

so that a firm in industry :math:`i` observes the noisy signals on that
:math:`\theta_t` presented to firms in both industries :math:`i` and
:math:`-i`.

The appropriate innovations representation becomes

.. math::
    :label: kf22&23

      \begin{aligned}
      \hat \theta_{t+1} & =   \rho
      \hat \theta_t + k a_t  \\
      w_t  & =  \begin{bmatrix} 1 \\ 1 \end{bmatrix} \hat \theta_t + a_t
     \end{aligned}

where :math:`a_t \equiv w_t - E [w_t | w^{t-1}]` is a
:math:`(2 \times 1)` vector of innovations in :math:`w_t` and :math:`k`
is now a :math:`(1 \times 2)` vector of Kalman gains.

Formulas for the Kalman filter imply that

.. math::
    :label: kf24

      \begin{aligned}
      k ={ \rho  p \over 2 p + \sigma_e^2}
      \begin{bmatrix}1 & 1 \end{bmatrix}
      \end{aligned}

where :math:`p = E \tilde \theta_t \tilde \theta_t^T` now satisfies the
Riccati equation

.. math::
    :label: ricc2

      \begin{aligned}
      p = \sigma_v^2 + {p \rho^2 \sigma_e^2 \over 2 p + \sigma_e^2}.
      \end{aligned}

Thus, when the representative firm in industry :math:`i` observes *two*
noisy signals on :math:`\theta_t`, we can express the equilibrium law of
motion for capital recursively as

.. math::
      :label: sol3

      \begin{aligned}
      k_{t+1}^i & =  \tilde \lambda k_t^i + {1 \over \lambda - \rho}\hat \theta_{t+1}  \\
      \hat \theta_{t+1} & =  \rho \theta_t + {\rho p \over 2  p + \sigma_e^2} (e_{1t}+e_{2t}) - {\rho \sigma_e^2 \over 2 p + \sigma_e^2} \tilde \theta_t \\
      \tilde \theta_{t+1} & =  { \rho \sigma_e^2 \over 2 p + \sigma_e^2} \tilde  \theta_t - {p \rho \over 2 p + \sigma_e^2}(e_{1t}+e_{2t}) +v_t \\
      \theta_{t+1} & =  \rho \theta_t + v_t . 
       \end{aligned}
:cite:`Pearlman_Sargent2005` verify that the above representation is equivalent
with what one obtains by using the machinery of :cite:`PCL`.

TOM REWRITE THE FOLLOWING
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We shall encounter versions of precisely these formulae again in section
:eq:`PCL2` where we compute the equilibrium of Townsend’s
model in which the representative firm in industry :math:`i` receives a
second noisy signal on :math:`\theta_t` by inferring it from
:math:`P_t^{-i}` and the other information that it has at time
:math:`t`. By extracting signals from the endogenous state variables, it
will turn out that the firm recovers exactly the same process for the
key additional state variable, the state reconstruction error
:math:`\tilde \theta_t`, that imperfect information contributes to the
dynamics.

System Description
==================


.. math::

   \begin{aligned}
      k_{t+1}^{i} & = & \tilde{\lambda}k_{t}^{i}+\frac{1}{\lambda-\rho}\hat{\theta}_{t+1}\\
      \hat{\theta}_{t+1} & = & \rho\theta_{t}+\frac{\rho p}{2p+\sigma_{e}^{2}}\left(e_{1,t}+e_{2,t}\right)-\frac{\rho\sigma_{e}^{2}}{2p+\sigma_{e}^{2}}\tilde{\theta}_{t}\\
      \tilde{\theta}_{t+1} & = & \frac{\rho\sigma_{e}^{2}}{2p+\sigma_{e}^{2}}\tilde{\theta}_{t}-\frac{p\rho}{2p+\sigma_{e}^{2}}\left(e_{1,t}+e_{2,t}\right)+v_{t}\\
      \theta_{t+1} & = & \rho\theta_{t}+v_{t}\\
      e_{1,t},e_{2,t} & \sim & \mathcal{N}\left(0,\sigma_{e}^{2}\right)\\
      v_{t} & \sim & \mathcal{N}\left(0,\sigma_{v}^{2}\right)
   \end{aligned}

where:


.. math::

   \begin{aligned}
      \left(\tilde{\lambda}-1\right)\left(\tilde{\lambda}-\frac{1}{\beta}\right) & = & b\tilde{\lambda}\\
      \left(\lambda-1\right)\left(\lambda-\frac{1}{\beta}\right) & = & b\lambda\\
      \tilde{\lambda} & \leq & \lambda\\
      p & = & \sigma_{v}^{2}+\frac{p\rho^{2}\sigma_{e}^{2}}{2p+\sigma_{e}^{2}}
   \end{aligned}

Parameters: :math:`\beta`, :math:`\rho`, :math:`b`, :math:`\sigma_v`,
and :math:`\sigma_e`

Computational Strategy
======================

Step 1: Solve for :math:`\tilde{\lambda}` and :math:`\lambda`
----------------------------------------------------------------

1. Cast
   :math:`\left(\lambda-1\right)\left(\lambda-\frac{1}{\beta}\right)=b\lambda`
   as :math:`p\left(\lambda\right)=0` where :math:`p` is a polynomial
   function of :math:`\lambda`.
2. Use ``numpy.roots`` to solve for the roots of :math:`p`
3. Verify :math:`\lambda \approx \frac{1}{\beta\tilde{\lambda}}`

Note that
:math:`p\left(\lambda\right)=\lambda^{2}-\left(1+b+\frac{1}{\beta}\right)\lambda+\frac{1}{\beta}`.

Step 2: Solve for :math:`p`
------------------------------

1. Cast
   :math:`p=\sigma_{v}^{2}+\frac{p\rho^{2}\sigma_{e}^{2}}{2p+\sigma_{e}^{2}}`
   as a discrete matrix Riccati equation.
2. Use ``quantecon.solve_discrete_riccati`` to solve for :math:`p`
3. Verify
   :math:`p \approx\sigma_{v}^{2}+\frac{p\rho^{2}\sigma_{e}^{2}}{2p+\sigma_{e}^{2}}`

Note that:


.. math::

   \begin{aligned}
      A & = & \left[\begin{array}{c}
      \rho\end{array}\right]\\
      B & = & \left[\begin{array}{c}
      \sqrt{2}\end{array}\right]\\
      R & = & \left[\begin{array}{c}
      \sigma_{e}^{2}\end{array}\right]\\
      Q & = & \left[\begin{array}{c}
      \sigma_{v}^{2}\end{array}\right]\\
      N & = & \left[\begin{array}{c}
      0\end{array}\right]
   \end{aligned}

Step 3: Represent the system using ``quantecon.LinearStateSpace``
-------------------------------------------------------------------

.. math::


   \begin{aligned}
      \left[\begin{array}{c}
      k_{t+1}^{i}\\
      \hat{\theta}_{t+1}\\
      \tilde{\theta}_{t+1}\\
      \theta_{t+1}
      \end{array}\right] & = & \underbrace{\left[\begin{array}{cccc}
      \tilde{\lambda} & 0 & \frac{1}{\lambda-\rho}\frac{-\rho\sigma_{e}^{2}}{2p+\sigma_{e}^{2}} & \frac{\rho}{\lambda-\rho}\\
      0 & 0 & \frac{-\rho\sigma_{e}^{2}}{2p+\sigma_{e}^{2}} & \rho\\
      0 & 0 & \frac{\rho\sigma_{e}^{2}}{2p+\sigma_{e}^{2}} & 0\\
      0 & 0 & 0 & \rho
      \end{array}\right]}_{A}\left[\begin{array}{c}
      k_{t}^{i}\\
      \hat{\theta}_{t}\\
      \tilde{\theta}_{t}\\
      \theta_{t}
      \end{array}\right]+\underbrace{\left[\begin{array}{ccc}
      \frac{\sigma_{e}}{\lambda-\rho}\frac{\rho p}{2p+\sigma_{e}^{2}} & \frac{\sigma_{e}}{\lambda-\rho}\frac{\rho p}{2p+\sigma_{e}^{2}} & 0\\
      \sigma_{e}\frac{\rho p}{2p+\sigma_{e}^{2}} & \sigma_{e}\frac{\rho p}{2p+\sigma_{e}^{2}} & 0\\
      -\sigma_{e}\frac{\rho p}{2p+\sigma_{e}^{2}} & -\sigma_{e}\frac{\rho p}{2p+\sigma_{e}^{2}} & \sigma_{v}\\
      0 & 0 & \sigma_{v}
      \end{array}\right]}_{C}\left[\begin{array}{c}
      z_{1,t+1}\\
      z_{2,t+1}\\
      z_{3,t+1}
      \end{array}\right]\\
      G & = & \left[\begin{array}{cccc}
      0 & 0 & 0 & 0\end{array}\right]\\
      H & = & \left[\begin{array}{c}
      0\end{array}\right]\\
      \left[\begin{array}{c}
      z_{1,t+1}\\
      z_{2,t+1}\\
      z_{3,t+1}
      \end{array}\right] & \sim & \mathcal{N}\left(0,I\right)
   \end{aligned}

Initial state:
:math:`\left[\begin{array}{ccccc} 0 & 0 & 0 & 0 \end{array}\right]'`

As usual, this representation is one of many possible representations.

System 1
========

.. math::

   \begin{aligned}
      \left[\begin{array}{c}
      k_{t+1}^{i}\\
      e_{t}\\
      \tilde{\theta}_{t+1}\\
      \theta_{t+1}
      \end{array}\right] & = & \underbrace{\left[\begin{array}{cccc}
      \tilde{\lambda} & 0 & \frac{1}{\lambda-\rho}\frac{-\rho\sigma_{e}^{2}}{2p+\sigma_{e}^{2}} & \frac{\rho}{\lambda-\rho}\\
      0 & 0 & 0 & 0\\
      0 & 0 & \frac{\rho\sigma_{e}^{2}}{2p+\sigma_{e}^{2}} & 0\\
      0 & 0 & 0 & \rho
      \end{array}\right]}_{A}\left[\begin{array}{c}
      k_{t}^{i}\\
      e_{t-1}\\
      \tilde{\theta}_{t}\\
      \theta_{t}
      \end{array}\right]+\underbrace{\left[\begin{array}{cc}
      \frac{\sigma_{e}}{\lambda-\rho}\frac{\rho p}{2p+\sigma_{e}^{2}} & 0\\
      \sigma_{e} & 0\\
      -\sigma_{e}\frac{\rho p}{2p+\sigma_{e}^{2}} & \sigma_{v}\\
      0 & \sigma_{v}
      \end{array}\right]}_{C}\left[\begin{array}{c}
      z_{1,t+1}\\
      z_{2,t+1}
      \end{array}\right]\\
      G & = & \left[\begin{array}{cccc}
      b & 1 & 0 & 1\end{array}\right]\\
      H & = & \left[\begin{array}{c}
      0\end{array}\right]\\
      \left[\begin{array}{c}
      z_{1,t+1}\\
      z_{2,t+1}
      \end{array}\right] & \sim & \mathcal{N}\left(0,I\right)
   \end{aligned}

Initial state:
:math:`\left[\begin{array}{ccccc} 0 & 0 & 0 & 0 \end{array}\right]'`

.. code-block:: ipython

    import numpy as np
    import quantecon as qe
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    import plotly.express as px
    import plotly.offline as pyo
    
    
    pyo.init_notebook_mode(connected=True)

.. code-block:: python3

    β = 0.9  # Discount factor
    ρ = 0.8  # Persistence parameter for the hidden state
    b = 0.5  # Demand curve parameter
    σ_v = 0.5  # Standard deviation of shock to θ_t 
    σ_e = 0.6  # Standard deviation of shocks to w_t

.. code-block:: python3

    # Compute λ
    poly = np.array([1, -(1 + β + b) / β, 1 / β])
    roots_poly = np.roots(poly)
    λ_tilde = roots_poly.min()
    λ = roots_poly.max()

.. code-block:: python3

    # Verify that λ = (βλ_tilde) ^ (-1)
    tol = 1e-12
    np.max(np.abs(λ - 1 / (β * λ_tilde))) < tol

.. code-block:: python3

    A_ricc = np.array([[ρ]])
    B_ricc = np.array([[np.sqrt(2)]])
    R_ricc = np.array([[σ_e ** 2]])
    Q_ricc = np.array([[σ_v ** 2]])
    N_ricc = np.zeros((1, 1))
    p = qe.solve_discrete_riccati(A_ricc, B_ricc, Q_ricc, R_ricc, N_ricc).item()

.. code-block:: python3

    # Verify that p = σ_v^2 + (pρ^2σ_e^2) / (2p + σ_e^2)
    tol = 1e-12
    np.abs(p - (σ_v ** 2 + p * ρ ** 2 * σ_e ** 2 / (2 * p + σ_e ** 2))) < tol

.. code-block:: python3

    term_0 = -ρ * σ_e ** 2 / (2 * p + σ_e ** 2)
    term_1 = ρ * p / (2 * p + σ_e ** 2)
    
    A_lss = np.array([[λ_tilde, 0., term_0 / (λ - ρ), ρ / (λ - ρ)],
                     [0., 0., 0., 0.],
                     [0., 0., -term_0, 0.,],
                     [0., 0., 0., ρ],])
    
    C_lss = np.array([[term_1 * σ_e / (λ - ρ),  0.],
                     [σ_e, 0.],
                     [-term_1 * σ_e,  σ_v],
                     [0., σ_v]])
    
    G_lss = np.array([[b, 1., 0., 1.]])

.. code-block:: python3

    mu_0 = np.array([0., 0., 0., 0.])
    
    lss = qe.LinearStateSpace(A_lss, C_lss, G_lss, mu_0=mu_0)

.. code-block:: python3

    ts_length = 100_000
    x, y = lss.simulate(ts_length)

.. code-block:: python3

    # Plot sample time path
    t = 300
    
    subplot_titles = [r'$k^{i}_t$',
                      r'$e_t$',
                      r'$\tilde{\theta}_t$',
                      r'$\theta_t$']
    
    fig = make_subplots(rows=x.shape[0], cols=1, subplot_titles=subplot_titles)
    
    for idx in range(x.shape[0]):
        fig.add_trace(go.Scatter(y=x[idx, :t],
                             legendgroup='trend'),
                 row=idx+1,
                 col=1)
        
    fig.update_layout(height=1200)    
        
    fig.show()

.. code-block:: python3

    px.line(y[0, :t])

.. code-block:: python3

    fig = px.histogram(x[2], title=r'$\mathrm{Histogram: }\: \tilde{\theta}_{t}$')
    fig.update_layout(height=500)   

.. code-block:: python3

    # Compute the mean of \tilde{\theta}
    x[2].mean()

System 2
========

.. math::

   \begin{aligned}
      \left[\begin{array}{c}
      k_{t+1}^{i}\\
      e_{1,t}\\
      e_{2,t}\\
      \tilde{\theta}_{t+1}\\
      \theta_{t+1}
      \end{array}\right] & = & \underbrace{\left[\begin{array}{ccccc}
      \tilde{\lambda} & 0 & 0 & \frac{1}{\lambda-\rho}\frac{-\rho\sigma_{e}^{2}}{2p+\sigma_{e}^{2}} & \frac{\rho}{\lambda-\rho}\\
      0 & 0 & 0 & 0 & 0\\
      0 & 0 & 0 & 0 & 0\\
      0 & 0 & 0 & \frac{\rho\sigma_{e}^{2}}{2p+\sigma_{e}^{2}} & 0\\
      0 & 0 & 0 & 0 & \rho
      \end{array}\right]}_{A}\left[\begin{array}{c}
      k_{t}^{i}\\
      e_{1,t-1}\\
      e_{2,t-1}\\
      \tilde{\theta}_{t}\\
      \theta_{t}
      \end{array}\right]+\underbrace{\left[\begin{array}{ccc}
      \frac{\sigma_{e}}{\lambda-\rho}\frac{\rho p}{2p+\sigma_{e}^{2}} & \frac{\sigma_{e}}{\lambda-\rho}\frac{\rho p}{2p+\sigma_{e}^{2}} & 0\\
      \sigma_{e} & 0 & 0\\
      0 & \sigma_{e} & 0\\
      -\sigma_{e}\frac{\rho p}{2p+\sigma_{e}^{2}} & -\sigma_{e}\frac{\rho p}{2p+\sigma_{e}^{2}} & \sigma_{v}\\
      0 & 0 & \sigma_{v}
      \end{array}\right]}_{C}\left[\begin{array}{c}
      z_{1,t+1}\\
      z_{2,t+1}\\
      z_{3,t+1}
      \end{array}\right]\\
      G & = & \left[\begin{array}{ccccc}
      b & 1 & 0 & 1 & 0\\
      b & 1 & 0 & 0 & 1
      \end{array}\right]\\
      H & = & \left[\begin{array}{c}
      0\\
      0
      \end{array}\right]\\
      \left[\begin{array}{c}
      z_{1,t+1}\\
      z_{2,t+1}\\
      z_{3,t+1}
      \end{array}\right] & \sim & \mathcal{N}\left(0,I\right)
   \end{aligned}

.. code-block:: python3

    term_0 = -ρ * σ_e ** 2 / (2 * p + σ_e ** 2)
    term_1 = ρ * p / (2 * p + σ_e ** 2)
    
    A_lss = np.array([[λ_tilde, 0., 0., term_0 / (λ - ρ), ρ / (λ - ρ)],
                     [0., 0., 0., 0., 0.],
                     [0., 0., 0., -term_0, 0.,],
                     [0., 0., 0., 0., ρ],])
    
    C_lss = np.array([[term_1 * σ_e / (λ - ρ), term_1 * σ_e / (λ - ρ), 0.],
                     [σ_e, 0., 0.],
                    [0., σ_e, 0.],
                     [-term_1 * σ_e, -term_1 * σ_e, σ_v],
                     [0., 0., σ_v]])
    
    G_lss = np.array([[b, 1., 0., 1. ,0.],
                      [b, 1., 0., 0., 1.]])

Various inserts
---------------

Townsend wanted to assume that at time :math:`t` firms in industry
:math:`i` observe :math:`k_t^i, Y_t^i, P_t^i, (P^{-i})^t`, where
:math:`(P^{-i})^t` is the history of prices in the other market up to
time :math:`t`.

(Because that turned out to be too challenging, Townsend made an
alternative assumption that eased his calculations: that after a large
number :math:`S` of periods, firms in industry :math:`i` observe the
noisy signals received by firms in industry :math:`-i` :math:`S` periods
ago and earlier. **TOM: ACTUALLY, I THINK HE ASSUMES THAT THEY SEE**
:math:`\theta_{t-S}`)

Request 1 for Quentin: Dec 12 2020
----------------------------------

Hi. Please take the above system and do the following things with it.

-  write it in state-space form with the state being
   :math:`k_t^i, \theta_t, \tilde \theta_t` with shock vector having
   components :math:`e_t, v_t`

-  add a *measurement equation* for
   :math:`P_t^i = b k_t^i + \theta_t + e_t` and :math:`\theta_t + e_t`
   and :math:`e_t`.

-  note that this is a system in which the state vector transition
   equation and the measurement equation have a common component. This
   is fine with the :math:`A,B,C,D` notation that LPH and I often use. I
   can explain details on zoom if you wish.

-  use a quantecon linear state space program to get impulse response
   functions for :math:`k_t^i` with respect to shocks :math:`v_t, e_t`.

-  use RMT5 chapter 2 page 45 or so formulas to compute covariance and
   cross-covariance matrices for the state and measurement vectors.

-  use formulas for multivariate normal distribution from RMT5 chapter 2
   to compute the population regression (and :math:`R^2`) of :math:`e_t`
   against :math:`P_t^i, k_t^i, \tilde \theta_t`

**note:** Maybe we need a helper function to compute the chapter 2
formulas for second moments for a state space system :math:`A, B, C, D`
– I can sketch the formulas for you. Zejin will be interested in these
formulas too. We’ll use them in *request 2* below too.

end of Quentin request 1
------------------------

Request 2 for Quentin: Dec 12 2020
----------------------------------

This request is an analogue of request 1 except for the two-noisy signal
model. Please take the above system and do the following things with it.

-  write it in state-space form with the state being
   :math:`k_t^i, \theta_t, \tilde \theta_t` with shock vector having
   components :math:`e_{1t}, e_{2t}, v_t`

-  add a *measurement equation* for
   :math:`P_t^i = b k_t^i + \theta_t + e_{it}` and
   :math:`\theta_t + e_{it}` and :math:`e_{it}` for our two
   :math:`i`\ s, :math:`i=1,2`.

-  note that this is a system in which the state vector transition
   equation and the measurement equation have common components. This is
   fine with the :math:`A,B,C,D` notation that LPH and I often use. I
   can explain details on zoom if you wish.

-  use a quantecon linear state space program to get impulse response
   functions for :math:`k_t^i` with respect to shocks
   :math:`v_t, e_{1t} + e_{2t}`.

-  use RMT5 chapter 2 page 45 or so formulas to compute covariance and
   cross-covariance matrices for the state and measurement vectors.

-  use formulas for multivariate normal distribution from RMT5 chapter 2
   to compute the population regression (and :math:`R^2`) of
   :math:`e_{-i,t}` against :math:`P_t^i, k_t^i, \tilde \theta_t` say
   for :math:`i=1`.

end of Quentin request 2
------------------------



.. rubric:: Footnotes

.. [#footnote1]  See :cite:`AHMS1996` for an account of invariant subspace methods. :cite:`ahms`

.. [#footnote2]  See :cite:`AHMS1996` for a discussion
      of the information assumptions needed to create a situation
      in which higher order beliefs appear in equilibrium decision rules.  The way
      to read our findings in light of :cite:`ams` is that Townsend's
      section 8 model  has too few sources of random shocks relative
      to sources of signals to permit higher order beliefs to
      play a role.

.. [#footnote3]  See :cite:`Sargent1987`, especially
      chapters IX and XIV, for the methods used in this section.

.. [#footnote4]  As noted :cite:`Sargent1987`, this difference equation is the Euler equation for
      the planning problem   of maximizing the discounted sum of consumer plus
      producer surplus.