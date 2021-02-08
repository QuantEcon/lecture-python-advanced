.. _knowing_the_forecast_of_others_v3:

.. include:: /_static/includes/header.raw

.. highlight:: python3


**************************************
Knowing the Forecasts of Others
**************************************


.. contents:: :depth: 2

In addition to what's in Anaconda, this lecture will need the following libraries:

.. code-block:: ipython
    :class: hide-output

    !pip install --upgrade quantecon
    !conda install -y -c plotly plotly plotly-orca

Introduction
=============

Robert E. Lucas, Jr. :cite:`lucas75`, Kenneth Kasa  :cite:`kasa`, and Robert Townsend
:cite:`townsend` showed that giving decision makers  incentives to infer persistent hidden  state
variables from equilibrium prices and quantities can
elongate and amplify impulse responses to aggregate 
shocks in business cycle models.

Townsend :cite:`townsend`
noted  that  such incentives can naturally
induce decision makers to want  to forecast the forecast of others.

This theme has been pursued and extended in analyses in which
decision makers' imperfect information forces them into pursuing an
infinite recursion of forming beliefs about the beliefs of other
(e.g., :cite:`ams`).

Lucas :cite:`lucas75` side stepped having decision makers forecast the
forecasts of other decision makers by assuming that they simply pool their
information before forecasting.

A **pooling equilibrium** like Lucas's plays a prominent role in this lecture.

Because he didn’t assume such pooling, :cite:`townsend`
confronted the forecasting the forecasts of others problem.

To formulate the  problem recursively required that Townsend define  decision maker's **state** vector.

Townsend concluded that his original  model required an intractable infinite dimensional  state space. 

Therefore, he constructed a more manageable approximating model in which the hidden Markov component of 
the demand shock is
revealed to all firms after a fixed and finite number of periods. 

In this lecture, as yet another instance of the theme that **finding the state is an art**, 
we show how to formulate Townsend's original model in terms of a low-dimensional state space.

By doing so, we show that Townsend's model shares equilibrium prices and quantities with those   that 
prevail in a pooling equilibrium.  

That finding emerged from a line of research about  Townsend's model that  culminated in  
:cite:`Pearlman_Sargent2005` that built on :cite:`PCL`.

However, rather than deploying the :cite:`PCL` machinery here, we shall rely instead on a sneaky
**guess-and-verify** tactic.

* We compute a pooling equilibrium and represent it as an instance of  a linear  state-space system provided by 
  the Python class ``quantecon.LinearStateSpace``.

* Leaving the state-transition equation for the pooling equilibrium unaltered, we alter the observation vector 
  for a firm to what it in in Townsend's original model. So rather than directly observing the signal received by 
  firms in the other industry, a firm sees the equilibrium price
  of the good produced by the other industry.

* We compute a population linear least squares regression of the noisy signal that firms in the other 
  industry receive in a pooling equilibrium on time :math:`t` information that a firm receives in Townsend's 
  original model. The :math:`R^2` in this regression equals :math:`1`.  That verifies that a firm's information 
  set in Townsend's original model equals its information set in a pooling equilibrium. Therefore, equilibrium
  prices and quantities in Townsend's original model equal those in a pooling equilibrium. 


A Sequence of Models
--------------------

We proceed by describing  a sequence of models of two industries that are linked in a
single way: shocks to the demand curves for their products have a common
component.

The models are simplified versions   of Townsend's  :cite:`townsend`.

Townsend's is a model of a rational expectations equilibrium in which firms confront
the problem  **forecasting the forecasts of others**.

In Townsend's model, firms condition  their forecasts on observed endogenous variables whose equilibrium laws of motion
are determined by their own forecasting functions.

We start with model  components that we shall progressively assemble in  ways that can help us to appreciate the structure of a
**pooling equilibrium**  that ultimately concerns us.

While keeping other aspects of the model the same, we shall study
consequences of alternative assumptions about what decision makers
observe.

Technically, this lecture deploys concepts and tools that appear
in  `First Look at Kalman Filter <https://python-intro.quantecon.org/kalman.html>`__ and
`Rational Expectations Equilibrium <https://python-intro.quantecon.org/rational_expectations.html>`__.



The Setting
============

We cast all variables in terms of deviations from means.

Therefore,  we omit constants from inverse demand curves
and other functions.

Firms in each of two industries :math:`i=1,2` use a single factor of
production, capital :math:`k_t^i`, to produce output of a single good,
:math:`y_t^i`.

Firms bear quadratic costs of adjusting their capital stocks.

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

The presence of costs of adjusting their capital stocks imparts to firms an incentives to forecast the price of the
good that they sell.

Throughout, we use the **rational expectations** equilibrium concept presented
in this lecture `Rational Expectations Equilibrium <https://python-intro.quantecon.org/rational_expectations.html>`__.

We let capital letters denote market wide objects and lower case letters
denote objects chosen by a representative firm.

In each industry, a competitive equilibrium prevails.

To rationalize the big :math:`K`, little :math:`k` connection, we can
think of there being a continua of each type of firm, each indexed by
:math:`\omega \in [0,1]` with
:math:`K^i = \int_0^1 k^i(\omega) d \omega`.

In equilibrium, :math:`k_t^i = K_t^i`, but as usual we must distinguish
between :math:`k_t^i` and :math:`K_t^i` when we pose the firm’s
optimization problem.


Tactics
===========

We shall  compute
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
   :math:`\{\theta_t, \epsilon_{t}^i\}` separately are never observed. However, at
   time :math:`t`, a history :math:`w^t` of a scalar noise-ridden
   observations on :math:`\theta_t` is observed at time :math:`t`.

-  **Two noise-ridden observations on** :math:`\theta_t`: Values of
   :math:`\{\theta_t, \epsilon_{t}^i\}` separately are never observed. However, at
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
:math:`t` but future :math:`\theta_{t+j}` and
:math:`\epsilon_{t+j}^i` are not observed.

To find an equilibrium when only a history :math:`w_t` of a single noise
ridden observations on :math:`\theta_t` is observed, we again apply a
certainty equivalence principle and replace future values of the random
variables :math:`\theta_s, \epsilon_{s}^i, s \geq t` with their
mathematical expectations conditioned on :math:`w^t`.

To find an equilibrium when only a history :math:`w_t` of a *two* noisy
signal on :math:`\theta_t` is observed, we replace future values of the
random variables :math:`\theta_s, \epsilon_{s}^i, s \geq t` with their
mathematical expectations conditioned on history :math:`w^t`.

We call the equilibrium with two noise-ridden observations on :math:`\theta_t` a **pooling equilibrium**.

-  It corresponds to an arrangement in which at the beginning of each
   period firms in industries :math:`1` and :math:`2` somehow get
   together and share information about current values of their noisy
   signals on :math:`\theta`.

We want ultimately to compare outcomes in  a pooling equilibrium
with an equilibrium under the following alternative information structure for a firm
in industry :math:`i` that interested :cite:`townsend`:

-  **Firm** :math:`i`\ ’s **noise-ridden signal on** :math:`\theta_t` **and the
   price in industry** :math:`-i`**:** At time :math:`t`, a firm in industry
   :math:`i` observes a history :math:`w^t` of *one* noise-ridden signal
   on :math:`\theta_t` and a history of industry :math:`-i`\ ’s price is
   observed.


With this information structure,  the representative firm :math:`i` sees the price as well as the
aggregate state variable :math:`Y_t^i` in its own industry.

That allows it to  infer
the total demand shock :math:`\theta_t + \epsilon_{t}^i`.

However, at time :math:`t`, the firm sees only :math:`P_t^{-i}` and does
not see :math:`Y_t^{-i}`, so that firm :math:`i` does not directly
observe :math:`\theta_t + \epsilon_t^{-i}`.

Nevertheless, it will turn out that   equilibrium prices and quantities in this equilibrium equal
their counterparts in a pooling equilibrium because
firms in industry :math:`i` are able to infer the noisy signal about the demand shock
received by firms in industry :math:`-i`.

We shall eventually verify this assertion  by using a guess and verify tactic. [#footnote0]_


Equilibrium conditions
=======================

It is convenient to solve the firm’s problem without
uncertainty  by forming the Lagrangian:

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

We can write  Euler equation :eq:`town4` as:


.. math::
    :label: pcl10

      \begin{aligned}
      g_t^i = P_t^i + \beta g_{t+1}^i
      \end{aligned}

In addition, we have the law of motion for :math:`\theta_t`,
:eq:`town2a`, and the demand equation
:eq:`town2`.

In summary, with perfect foresight, equilibrium conditions for industry
:math:`i` include the following system of difference equations:




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

Equilibrium under perfect foresight
------------------------------------

Our first step is to compute the equilibrium law of motion for
:math:`k_t^i` under perfect foresight.

Let :math:`L` be the lag
operator. [#footnote3]_

Equations :eq:`pcl10` and :eq:`pcl11`
imply the second order difference equation in
:math:`k_t^i`: [#footnote4]_


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
      (\epsilon_{t+1}^i + \theta_{t+1}  ) .
      \end{aligned}

Recall that we have already set :math:`k^i = K^i` at the appropriate point in the argument (i.e., _after_ having derived the first-order necessary
conditions for a representative firm in industry :math:`i`.

Thus,  under perfect foresight the equilibrium capital stock in industry :math:`i` satisfies

.. math::
    :label: town5

      \begin{aligned}
      k_{t+1}^i = \tilde \lambda k_t^i + \sum_{j=1}^\infty (\tilde \lambda \beta)^j
      (\epsilon_{t+j}^i +  \theta_{t+j}) .
      \end{aligned}

Next, we shall investigate consequences of replacing future  values of
:math:`(\epsilon_{t+j}^i +  \theta_{t+j})` in equation :eq:`town5`  with alternative forecasting schemes.

In particular, we shall compute   equilibrium laws of motion for capital
under alternative assumptions about the information available to
decision makers in market :math:`i`.

Equilibrium with :math:`\theta_t` stochastic but observed at :math:`t`
=======================================================================

If future :math:`\theta`\ ’s are unknown at :math:`t`, it is appropriate
to replace all random variables on the right side of
:eq:`town5` with their conditional expectations based on
the information available to decision makers in market :math:`i`.

For now, we assume that this information set
:math:`I_t^p = \begin{bmatrix} \theta^t & \epsilon^{it} \end{bmatrix}`,
where :math:`z^t` represents the infinite history of variable
:math:`z_s` up to time :math:`t`.

Later we shall give firms less information.

To obtain an appropriate counterpart to :eq:`town5` under our
current assumption about information, we apply a certainty equivalence
principle.

In particular, it is appropriate to take :eq:`town5` and
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

For future purposes, it is useful to represent the equilibrium
:math:`\{k_t^i\}_t` process recursively as


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

We get closer to a model that we ultimately want to study by now
assuming that firms in market :math:`i` do not observe :math:`\theta_t`,
but instead observe a history  :math:`w^t` of noisy signals at time :math:`t`.

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

where :math:`w^t = [w_t, w_{t-1}, \ldots, w_0]` denotes the history of the :math:`w_s` process up to
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


:math:`\theta`-reconstruction error:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Define the state *reconstruction error* :math:`\tilde \theta_t` by

.. math::

   \begin{aligned}
   \tilde \theta_t = \theta_t - \hat \theta_t .
   \end{aligned}

Then :math:`p = E \tilde \theta_t^2`.

Equations :eq:`kf1&2` and :eq:`kf3&4` imply

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
the equilibrium when :math:`w^t` is observed.

Relative to :eq:`solution1`, the equilibrium acquires a new
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
      :label: sol0a

      \begin{aligned}
      k_{t+1}^i & =  \tilde \lambda k_t^i + {1 \over \lambda - \rho}   \hat \theta_{t+1}  \\
      \hat \theta_{t+1} & =  \rho \theta_t + {\rho p \over p + \sigma_e^2} e_t - {\rho \sigma_e^2 \over p + \sigma_e^2} \tilde \theta_t  \\
      \tilde \theta_{t+1} & = { \rho \sigma_e^2 \over p + \sigma_e^2} \tilde
         \theta_t - {p \rho \over p + \sigma_e^2} e_t + v_t  \\
      \theta_{t+1} & =  \rho \theta_t + v_t .  \end{aligned}

Two noisy signals
--------------------

We now construct a **pooling equilibrium** by assuming that a firm in
industry :math:`i` receives a vector :math:`w_t` of *two* noisy signals
on :math:`\theta_t`:


.. math::

      \begin{eqnarray*}
      \theta_{t+1} & = & \rho\theta_{t}+v_{t} \label{kf20} \\
      w_{t} & = & \begin{bmatrix}1\\
      1
      \end{bmatrix}\theta_{t}+\begin{bmatrix}e_{1t} \\
      e_{2t}
      \end{bmatrix} \label{kf21}
      \end{eqnarray*}

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

Thus, when a representative firm in industry :math:`i` observes *two*
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


Below, by using a guess-and-verify tactic,  we shall show that outcomes in this  **pooling equilibrium** equal those in an equilibrium under the alternative
information structure that interested :cite:`townsend`. [#footnote5]_


Guess-and-verify tactic
===================================

As a preliminary step we shall take our recursive representation :eq:`sol0a`
of an equilibrium in industry :math:`i` with one noisy signal
on :math:`\theta_t` and perform the following steps:

- Compute :math:`\lambda` and :math:`\tilde{\lambda}` by posing a
  root-finding problem and then solving it using ``numpy.roots``

- Compute :math:`p` by forming the appropriate  discrete Riccati equation and then solving it
  using ``quantecon.solve_discrete_riccati``

- Add a *measurement equation* for
  :math:`P_t^i = b k_t^i + \theta_t + e_t`, :math:`\theta_t + e_t`,
  and :math:`e_t` to system :eq:`sol0a`. Write the resulting system
  in state-space form and encode it using ``quantecon.LinearStateSpace``

- Use methods of the ``quantecon.LinearStateSpace`` to compute impulse response
  functions of :math:`k_t^i` with respect to shocks :math:`v_t, e_t`.

After analyzing the one-noisy-signal structure in this way,  by making appropriate modifications
we shall analyze the two-noisy-signal
structure.

We proceed to analyze first the one-noisy-signal structure and then the two-noisy-signal structure.


Equilibrium with one signal on :math:`\theta_t`
================================================

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

We use the following representation for constructing the
``quantecon.LinearStateSpace`` instance.

.. math::

    \begin{eqnarray*}
    \underbrace{\left[\begin{array}{c}
    e_{t+1}\\
    k_{t+1}^{i}\\
    \tilde{\theta}_{t+1}\\
    P_{t+1}\\
    \theta_{t+1}\\
    v_{t+1}
    \end{array}\right]}_{x_{t+1}} & = & \underbrace{\left[\begin{array}{cccccc}
    0 & 0 & 0 & 0 & 0 & 0\\
    \frac{\kappa}{\lambda-\rho} & \tilde{\lambda} & \frac{-1}{\lambda-\rho}\frac{\kappa\sigma_{e}^{2}}{p} & 0 & \frac{\rho}{\lambda-\rho} & 0\\
    -\kappa & 0 & \frac{\kappa\sigma_{e}^{2}}{p} & 0 & 0 & 1\\
    \frac{b\kappa}{\lambda-\rho} & b\tilde{\lambda} & \frac{-b}{\lambda-\rho}\frac{\kappa\sigma_{e}^{2}}{p} & 0 & \frac{b\rho}{\lambda-\rho}+\rho & 1\\
    0 & 0 & 0 & 0 & \rho & 1\\
    0 & 0 & 0 & 0 & 0 & 0
    \end{array}\right]}_{A}\underbrace{\left[\begin{array}{c}
    e_{t}\\
    k_{t}^{i}\\
    \tilde{\theta}_{t}\\
    P_{t}\\
    \theta_{t}\\
    v_{t}
    \end{array}\right]}_{x_{t}}+\underbrace{\left[\begin{array}{cc}
    \sigma_{e} & 0\\
    0 & 0\\
    0 & 0\\
    \sigma_{e} & 0\\
    0 & 0\\
    0 & \sigma_{v}
    \end{array}\right]}_{C}\left[\begin{array}{c}
    z_{1,t+1}\\
    z_{2,t+1}
    \end{array}\right]\\
    \underbrace{\left[\begin{array}{c}
    P_{t}\\
    e_{t}+\theta_{t}\\
    e_{t}
    \end{array}\right]}_{y_{t}} & = & \underbrace{\left[\begin{array}{cccccc}
    0 & 0 & 0 & 1 & 0 & 0\\
    1 & 0 & 0 & 0 & 1 & 0\\
    1 & 0 & 0 & 0 & 0 & 0
    \end{array}\right]}_{G}\underbrace{\left[\begin{array}{c}
    e_{t}\\
    k_{t}^{i}\\
    \tilde{\theta}_{t}\\
    P_{t}\\
    \theta_{t}\\
    v_{t}
    \end{array}\right]}_{x_{t}}+\underbrace{\left[\begin{array}{c}
    0\\
    0\\
    0
    \end{array}\right]}_{H}w_{t+1}\\
    \left[\begin{array}{c}
    z_{1,t+1}\\
    z_{2,t+1}\\
    w_{t+1}
    \end{array}\right] & \sim & \mathcal{N}\left(0,I\right)\\
    \kappa & = & \frac{\rho p}{p+\sigma_{e}^{2}}
    \end{eqnarray*}

This representation includes extraneous variables  such as :math:`P_{t}` in the
state vector.

We formulate things in this way because
it allows us easily to compute covariances of these variables with other
components of the state vector (step 5 above) by using the ``stationary_distributions`` method of the ``LinearStateSpace`` class.


.. code-block:: ipython

    import numpy as np
    import quantecon as qe
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    import plotly.express as px
    import plotly.offline as pyo
    from statsmodels.regression.linear_model import OLS
    from IPython.display import display, Latex, Image

    pyo.init_notebook_mode(connected=True)

.. code-block:: python3

    β = 0.9  # Discount factor
    ρ = 0.8  # Persistence parameter for the hidden state
    b = 1.5  # Demand curve parameter
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
    B_ricc = np.array([[1.]])
    R_ricc = np.array([[σ_e ** 2]])
    Q_ricc = np.array([[σ_v ** 2]])
    N_ricc = np.zeros((1, 1))
    p = qe.solve_discrete_riccati(A_ricc, B_ricc, Q_ricc, R_ricc, N_ricc).item()

    p_one = p  # Save for comparison later

.. code-block:: python3

    # Verify that p = σ_v ^ 2 + p * ρ ^ 2 - (ρ * p) ^ 2 / (p + σ_e ** 2)
    tol = 1e-12
    np.abs(p - (σ_v ** 2 + p * ρ ** 2 - (ρ * p) ** 2 / (p + σ_e ** 2))) < tol

.. code-block:: python3

    κ = ρ * p / (p + σ_e ** 2)
    κ_prod = κ * σ_e ** 2 / p

    κ_one = κ  # Save for comparison later

    A_lss = np.array([[0., 0., 0., 0., 0., 0.],
                     [κ / (λ - ρ), λ_tilde, -κ_prod / (λ - ρ), 0., ρ / (λ - ρ), 0.],
                     [-κ, 0., κ_prod, 0., 0., 1.],
                     [b * κ / (λ - ρ) , b * λ_tilde, -b * κ_prod / (λ - ρ), 0., b * ρ / (λ - ρ) + ρ, 1.],
                     [0., 0., 0., 0., ρ, 1.],
                     [0., 0., 0., 0., 0., 0.]])

    C_lss = np.array([[σ_e, 0.],
                     [0.,  0.],
                     [0.,  0.],
                     [σ_e, 0.],
                     [0., 0.],
                     [0., σ_v]])

    G_lss = np.array([[0., 0., 0., 1., 0., 0.],
                     [1., 0., 0., 0., 1., 0.],
                     [1., 0., 0., 0., 0., 0.]])

.. code-block:: python3

    mu_0 = np.array([0., 0., 0., 0., 0., 0.])

    lss = qe.LinearStateSpace(A_lss, C_lss, G_lss, mu_0=mu_0)

.. code-block:: python3

    ts_length = 100_000
    x, y = lss.simulate(ts_length, random_state=1)

.. code-block:: python3

    # Verify that two ways of computing P_t match
    np.max(np.abs(np.array([[1., b, 0., 0., 1., 0.]]) @ x - x[3])) < 1e-12

Step 4: Compute impulse response functions
-------------------------------------------------------------------

To compute impulse response functions of :math:`k_t^i`, we use the ``impulse_response`` method of the
``quantecon.LinearStateSpace`` class and plot the result.

.. code-block:: python3

    xcoef, ycoef = lss.impulse_response(j=21)
    data = np.array([xcoef])[0, :, 1, :]

    fig = go.Figure(data=go.Scatter(y=data[:-1, 0], name=r'$e_{t+1}$'))
    fig.add_trace(go.Scatter(y=data[1:, 1], name=r'$v_{t+1}$'))
    fig.update_layout(title=r'Impulse Response Function',
                       xaxis_title='Time',
                       yaxis_title=r'$k^{i}_{t}$')
    fig = fig1
    # Export to PNG file
    Image(fig1.to_image(format="png"))
    # fig1.show() will provide interactive plot when running
    # notebook locally

Step 5: Compute stationary covariance matrices and population regressions
-------------------------------------------------------------------------

We compute stationary covariance matrices  by
calling the ``stationary_distributions`` method of
the ``quantecon.LinearStateSpace`` class.

By appropriately decomposing the covariance matrix of the state vector, we obtain ingredients
of some population regression coefficients.

.. math::

  \Sigma_{x}=\left[\begin{array}{cc}
  \Sigma_{11} & \Sigma_{12}\\
  \Sigma_{21} & \Sigma_{22}
  \end{array}\right]

where :math:`\Sigma_{11}` is the covariance matrix of dependent variables and :math:`\Sigma_{22}` is the covariance matrix of independent variables.

Regression coefficients are  :math:`\beta=\Sigma_{21}\Sigma_{22}^{-1}`.

To verify an instance of a law of large numbers computation, we construct a long  simulation of 
the state vector and for the resulting sample compute
the ordinary least-squares estimator  of :math:`\beta` that we shall compare to the corresponding population regression
coefficients.

.. code-block:: python3

    _, _, Σ_x, Σ_y = lss.stationary_distributions()

    Σ_11 = Σ_x[0, 0]
    Σ_12 = Σ_x[0, 1:4]
    Σ_21 = Σ_x[1:4, 0]
    Σ_22 = Σ_x[1:4, 1:4]

    reg_coeffs = Σ_12 @ np.linalg.inv(Σ_22)

    print('Regression coefficients (e_t on k_t, P_t, \\tilde{\\theta_t})')
    print('------------------------------')
    print(r'k_t:', reg_coeffs[0])
    print(r'\tilde{\theta_t}:', reg_coeffs[1])
    print(r'P_t:', reg_coeffs[2])

.. code-block:: python3

    # Compute R squared
    R_squared = reg_coeffs @ Σ_x[1:4, 1:4] @ reg_coeffs  / Σ_x[0, 0]
    R_squared

.. code-block:: python3

    # Verify that the computed coefficients are close to least squares estimates
    model = OLS(x[0], x[1:4].T)
    reg_res = model.fit()
    np.max(np.abs(reg_coeffs - reg_res.params)) < 1e-2

.. code-block:: python3

    # Verify that R_squared matches least squares estimate
    np.abs(reg_res.rsquared - R_squared) < 1e-2

.. code-block:: python3

    # Verify that θ_t + e_t can be recovered
    model = OLS(y[1], x[1:4].T)
    reg_res = model.fit()
    np.abs(reg_res.rsquared - 1.) < 1e-6

Equilibrium with two noisy signals on :math:`\theta_t`
=========================================================

Steps 1, 4, and 5 are identical to those for the  one-noisy-signal structure.

Step 2  requires only a straightforward modification.

For step 3, we use construct the  following state-space representation so that we can get our hands on 
all of the random processes that we require in order  to compute a regression of the noisy signal about 
:math:`\theta` from the other industry that a firm receives directly in a pooling equilibrium on the information that
a firm receives in Townsend's original model.

For this purpose, we include  equilibrium goods prices from  both industries appear in the state vector:

.. math::

    \begin{eqnarray*}
    \underbrace{\left[\begin{array}{c}
    e_{1,t+1}\\
    e_{2,t+1}\\
    k_{t+1}^{i}\\
    \tilde{\theta}_{t+1}\\
    P_{t+1}^{1}\\
    P_{t+1}^{2}\\
    \theta_{t+1}\\
    v_{t+1}
    \end{array}\right]}_{x_{t+1}} & = & \underbrace{\left[\begin{array}{cccccccc}
    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
    \frac{\kappa}{\lambda-\rho} & \frac{\kappa}{\lambda-\rho} & \tilde{\lambda} & \frac{-1}{\lambda-\rho}\frac{\kappa\sigma_{e}^{2}}{p} & 0 & 0 & \frac{\rho}{\lambda-\rho} & 0\\
    -\kappa & -\kappa & 0 & \frac{\kappa\sigma_{e}^{2}}{p} & 0 & 0 & 0 & 1\\
    \frac{b\kappa}{\lambda-\rho} & \frac{b\kappa}{\lambda-\rho} & b\tilde{\lambda} & \frac{-b}{\lambda-\rho}\frac{\kappa\sigma_{e}^{2}}{p} & 0 & 0 & \frac{b\rho}{\lambda-\rho}+\rho & 1\\
    \frac{b\kappa}{\lambda-\rho} & \frac{b\kappa}{\lambda-\rho} & b\tilde{\lambda} & \frac{-b}{\lambda-\rho}\frac{\kappa\sigma_{e}^{2}}{p} & 0 & 0 & \frac{b\rho}{\lambda-\rho}+\rho & 1\\
    0 & 0 & 0 & 0 & 0 & 0 & \rho & 1\\
    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0
    \end{array}\right]}_{A}\underbrace{\left[\begin{array}{c}
    e_{1,t}\\
    e_{2,t}\\
    k_{t}^{i}\\
    \tilde{\theta}_{t}\\
    P_{t}^{1}\\
    P_{t}^{2}\\
    \theta_{t}\\
    v_{t}
    \end{array}\right]}_{x_{t}}+\underbrace{\left[\begin{array}{ccc}
    \sigma_{e} & 0 & 0\\
    0 & \sigma_{e} & 0\\
    0 & 0 & 0\\
    0 & 0 & 0\\
    \sigma_{e} & 0 & 0\\
    0 & \sigma_{e} & 0\\
    0 & 0 & 0\\
    0 & 0 & \sigma_{v}
    \end{array}\right]}_{C}\left[\begin{array}{c}
    z_{1,t+1}\\
    z_{2,t+1}\\
    z_{3,t+1}
    \end{array}\right]\\
    \underbrace{\left[\begin{array}{c}
    P_{t}^{1}\\
    P_{t}^{2}\\
    e_{1,t}+\theta_{t}\\
    e_{2,t}+\theta_{t}\\
    e_{1,t}\\
    e_{2,t}
    \end{array}\right]}_{y_{t}} & = & \underbrace{\left[\begin{array}{cccccccc}
    0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\
    0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\\
    1 & 0 & 0 & 0 & 0 & 0 & 1 & 0\\
    0 & 1 & 0 & 0 & 0 & 0 & 1 & 0\\
    1 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
    0 & 1 & 0 & 0 & 0 & 0 & 0 & 0
    \end{array}\right]}_{G}\underbrace{\left[\begin{array}{c}
    e_{1,t}\\
    e_{2,t}\\
    k_{t}^{i}\\
    \tilde{\theta}_{t}\\
    P_{t}^{1}\\
    P_{t}^{2}\\
    \theta_{t}\\
    v_{t}
    \end{array}\right]}_{x_{t}}+\underbrace{\left[\begin{array}{c}
    0\\
    0\\
    0\\
    0\\
    0\\
    0
    \end{array}\right]}_{H}w_{t+1}\\
    \left[\begin{array}{c}
    z_{1,t+1}\\
    z_{2,t+1}\\
    z_{3,t+1}\\
    w_{t+1}
    \end{array}\right] & \sim & \mathcal{N}\left(0,I\right)\\
    \kappa & = & \frac{\rho p}{2p+\sigma_{e}^{2}}
    \end{eqnarray*}

.. code-block:: python3

    A_ricc = np.array([[ρ]])
    B_ricc = np.array([[np.sqrt(2)]])
    R_ricc = np.array([[σ_e ** 2]])
    Q_ricc = np.array([[σ_v ** 2]])
    N_ricc = np.zeros((1, 1))
    p = qe.solve_discrete_riccati(A_ricc, B_ricc, Q_ricc, R_ricc, N_ricc).item()

    p_two = p  # Save for comparison later

.. code-block:: python3

    # Verify that p = σ_v^2 + (pρ^2σ_e^2) / (2p + σ_e^2)
    tol = 1e-12
    np.abs(p - (σ_v ** 2 + p * ρ ** 2 * σ_e ** 2 / (2 * p + σ_e ** 2))) < tol

.. code-block:: python3

    κ = ρ * p / (2 * p + σ_e ** 2)
    κ_prod = κ * σ_e ** 2 / p

    κ_two = κ  # Save for comparison later

    A_lss = np.array([[0., 0., 0., 0., 0., 0., 0., 0.],
                     [0., 0., 0., 0., 0., 0., 0., 0.],
                     [κ / (λ - ρ), κ / (λ - ρ), λ_tilde, -κ_prod / (λ - ρ), 0., 0., ρ / (λ - ρ), 0.],
                     [-κ, -κ, 0., κ_prod, 0., 0., 0., 1.],
                     [b * κ / (λ - ρ), b * κ / (λ - ρ), b * λ_tilde, -b * κ_prod / (λ - ρ), 0., 0., b * ρ / (λ - ρ) + ρ, 1.],
                     [b * κ / (λ - ρ), b * κ / (λ - ρ), b * λ_tilde, -b * κ_prod / (λ - ρ), 0., 0., b * ρ / (λ - ρ) + ρ, 1.],
                     [0., 0., 0., 0., 0., 0., ρ, 1.],
                     [0., 0., 0., 0., 0., 0., 0., 0.]])

    C_lss = np.array([[σ_e, 0., 0.],
                     [0., σ_e, 0.],
                     [0., 0.,  0.],
                     [0., 0.,  0.],
                     [σ_e, 0., 0.],
                     [0., σ_e, 0.],
                     [0., 0., 0.],
                     [0., 0., σ_v]])

    G_lss = np.array([[0., 0., 0., 0., 1., 0., 0., 0.],
                     [0., 0, 0, 0., 0., 1., 0., 0.],
                     [1., 0., 0., 0., 0., 0., 1., 0.],
                     [0., 1., 0., 0., 0., 0., 1., 0.],
                     [1., 0., 0., 0., 0., 0., 0., 0.],
                     [0., 1., 0., 0., 0., 0., 0., 0.]])


.. code-block:: python3

    mu_0 = np.array([0., 0., 0., 0., 0., 0., 0., 0.])

    lss = qe.LinearStateSpace(A_lss, C_lss, G_lss, mu_0=mu_0)

.. code-block:: python3

    ts_length = 100_000
    x, y = lss.simulate(ts_length, random_state=1)

.. code-block:: python3

    xcoef, ycoef = lss.impulse_response(j=20)

.. code-block:: python3

    data = np.array([xcoef])[0, :, 2, :]

    fig = go.Figure(data=go.Scatter(y=data[:-1, 0], name=r'$e_{1,t+1}$'))
    fig.add_trace(go.Scatter(y=data[:-1, 1], name=r'$e_{2,t+1}$'))
    fig.add_trace(go.Scatter(y=data[1:, 2], name=r'$v_{t+1}$'))
    fig.update_layout(title=r'Impulse Response Function',
                       xaxis_title='Time',
                       yaxis_title=r'$k^{i}_{t}$')
    fig2=fig
    # Export to PNG file
    Image(fig2.to_image(format="png"))
    # fig2.show() will provide interactive plot when running
    # notebook locally


.. code-block:: python3

    _, _, Σ_x, Σ_y = lss.stationary_distributions()

    Σ_11 = Σ_x[1, 1]
    Σ_12 = Σ_x[1, 2:5]
    Σ_21 = Σ_x[2:5, 1]
    Σ_22 = Σ_x[2:5, 2:5]

    reg_coeffs = Σ_12 @ np.linalg.inv(Σ_22)

    print('Regression coefficients (e_{2,t} on k_t, P^{1}_t, \\tilde{\\theta_t})')
    print('------------------------------')
    print(r'k_t:', reg_coeffs[0])
    print(r'\tilde{\theta_t}:', reg_coeffs[1])
    print(r'P_t:', reg_coeffs[2])

.. code-block:: python3

    # Compute R squared
    R_squared = reg_coeffs @ Σ_x[2:5, 2:5] @ reg_coeffs  / Σ_x[1, 1]
    R_squared

.. code-block:: python3

    # Verify that the computed coefficients are close to least squares estimates
    model = OLS(x[1], x[2:5].T)
    reg_res = model.fit()
    np.max(np.abs(reg_coeffs - reg_res.params)) < 1e-2

.. code-block:: python3

    # Verify that R_squared matches least squares estimate
    np.abs(reg_res.rsquared - R_squared) < 1e-2

.. code-block:: python3

    _, _, Σ_x, Σ_y = lss.stationary_distributions()

    Σ_11 = Σ_x[1, 1]
    Σ_12 = Σ_x[1, 2:6]
    Σ_21 = Σ_x[2:6, 1]
    Σ_22 = Σ_x[2:6, 2:6]

    reg_coeffs = Σ_12 @ np.linalg.inv(Σ_22)

    print('Regression coefficients (e_{2,t} on k_t, P^{1}_t, P^{2}_t, \\tilde{\\theta_t})')
    print('------------------------------')
    print(r'k_t:', reg_coeffs[0])
    print(r'\tilde{\theta_t}:', reg_coeffs[1])
    print(r'P^{1}_t:', reg_coeffs[2])
    print(r'P^{2}_t:', reg_coeffs[3])

.. code-block:: python3

    # Compute R squared
    R_squared = reg_coeffs @ Σ_x[2:6, 2:6] @ reg_coeffs  / Σ_x[1, 1]
    R_squared


## Key step

Now we come to the key step of verifying that equilibrium outcomes for prices and quantities are identical 
in the pooling equilibrium and Townsend's original model.

We accomplish this by  compute a population linear least squares regression of the noisy signal that firms in the other 
industry receive in a pooling equilibrium on time :math:`t` information that a firm receives in Townsend's 
original model.

Let's compute the regression and stare at the :math:`R^2`: 




.. code-block:: python3

    # Verify that θ_t + e^{2}_t can be recovered

    # θ_t + e^{2}_t on k^{i}_t, P^{1}_t, P^{2}_t, \\tilde{\\theta_t}
 

    model = OLS(y[1], x[2:6].T)
    reg_res = model.fit()
    np.abs(reg_res.rsquared - 1.) < 1e-6

.. code-block:: python3

    reg_res.rsquared



The :math:`R^2` in this regression equals :math:`1`. 

That verifies that a firm's information 
set in Townsend's original model equals its information set in a pooling equilibrium. 

Therefore, equilibrium prices and quantities in Townsend's original model equal those in a pooling equilibrium. 


Comparison of the two signal structures
======================================================================

It is enlightening side by side to  plot impulse response functions for capital in an industry for the two
information noisy-signal information structures.

Please remember that the two-signal structure corresponds to the **pooling equilibrium** and also
**Townsend's original model**.

.. code-block:: python3

    fig_comb = go.Figure(data=[*fig1.data,
                             *fig2.update_traces(xaxis='x2', yaxis='y2').data]).set_subplots(1, 2,
                                                                                             subplot_titles=("One noisy-signal structure", "Two noisy-signal structure"),
                                                                                             horizontal_spacing=0.1,
                                                                                             shared_yaxes=True)
    # Export to PNG file
    Image(fig_comb.to_image(format="png"))
    # fig_comb.show() will provide interactive plot when running
    # notebook locally


The graphs above show that

* the response of :math:`k_t^i` to  shocks :math:`v_t` to the hidden Markov demand state :math:`\theta_t` process  is **larger** in 
  the two-noisy=signal structure
* the response of :math:`k_t^i` to idiosyncratic *own-market*   noise-shocks :math:`e_t` is **smaller** in the two-noisy-signal structure


Taken together, these  findings in turn can be shown to imply that time series correlations and coherences between outputs in
the two industries are higher in the two-noisy-signals or **pooling** model.

The enhanced influence of the shocks :math:`v_t` to the hidden Markov demand state :math:`\theta_t` process that
emerges from the two-noisy-signal model relative to the one-noisy-signal model is a  symptom of a lower
equilibrium hidden-state  reconstruction error variance in the two-signal model:

.. code-block:: python3

    display(Latex('$\\textbf{Reconstruction error variances}$'))
    display(Latex(f'One-noise structure: ${round(p_one, 6)}$'))
    display(Latex(f'Two-noise structure: ${round(p_two, 6)}$'))

Kalman gains  for the two
structures are

.. code-block:: python3

    display(Latex('$\\textbf{Kalman Gains}$'))
    display(Latex(f'One noisy-signal structure: ${round(κ_one, 6)}$'))
    display(Latex(f'Two noisy-signals structure: ${round(κ_two, 6)}$'))


Notes on History of the Problem
=================================

To truncate what he saw as an intractable, infinite dimensional  state space,
Townsend constructed an approximating model in which the common hidden Markov demand shock
is revealed to all firms  after a fixed number of periods.

Thus,

  * Townsend wanted to assume that at time :math:`t` firms in industry :math:`i` observe
    :math:`k_t^i, Y_t^i, P_t^i, (P^{-i})^t`, where :math:`(P^{-i})^t` is the history of prices in
    the other market up to time :math:`t`.

  * Because that turned out to be too challenging, Townsend made an
    alternative assumption that eased his calculations: that after a large
    number :math:`S` of periods, firms in industry :math:`i` observe the
    hidden Markov component of the demand shock :math:`\theta_{t-S}`.

Townsend argued that the more manageable model  could do a good job of
approximating the intractable model in which the Markov component of the demand shock remains unobserved
for ever.

By applying technical machinery of :cite:`PCL`,
:cite:`Pearlman_Sargent2005` showed that there is a recursive
representation of the equilibrium of the perpetually and symmetrically
uninformed model formulated but not completely solved in section 8 of
:cite:`townsend`.

A reader of :cite:`Pearlman_Sargent2005` will notice that their representation of the equilibrium of 
Townsend's model exactly matches that of the  **pooling equilibrium** presented here.  

We have structured  our notation in  this lecture to faciliate comparison of the **pooling equilibrium**
constructed here with the equilibrium of Townsend's model reported in  :cite:`Pearlman_Sargent2005`.

The computational method of :cite:`Pearlman_Sargent2005` is recursive:
it enlists the Kalman filter and invariant subspace methods for
solving systems of Euler
equations [#footnote1]_ .

As :cite:`singleton`,
:cite:`kasa`, and :cite:`sargent91` also
found, the equilibrium is fully revealing: observed prices tell
participants in industry :math:`i` all of the information held by
participants in market :math:`-i` (:math:`-i` means not :math:`i`).

This
means that higher-order beliefs play no role: seeing equilibrium prices
in effect lets decision makers pool their information
sets [#footnote2]_ .

The disappearance of higher order beliefs means that
decision makers in this model do not really face a problem of
forecasting the forecasts of others. 

They know those forecasts because
they are the same as their own.

Further historical remarks
--------------------------

:cite:`sargent91` proposed a way to compute an equilibrium
without making Townsend’s approximation.

Extending the reasoning of :cite:`muth1960`, Sargent noticed that it is possible to
summarize the relevant history with a low dimensional object, namely, a
small number of current and lagged forecasting errors.

Positing an
equilibrium in a space of perceived laws of motion for endogenous
variables that takes the form of a vector autoregressive, moving
average, Sargent described an equilibrium as a fixed point of a mapping
from the perceived law of motion to the actual law of motion of that
form.

Sargent worked in the time domain and had to guess and verify the
appropriate orders of the autoregressive and moving average pieces of
the equilibrium representation.

By working in the frequency
domain :cite:`kasa` showed how to discover the appropriate
orders of the autoregressive and moving average parts, and also how to
compute an equilibrium.

The  :cite:`Pearlman_Sargent2005` recursive computational method, which stays in the time domain, also
discovered appropriate orders of the autoregressive and moving
average pieces.

In addition, by displaying equilibrium representations
in the form of :cite:`PCL`, :cite:`Pearlman_Sargent2005`
showed how the moving average piece is linked to the innovation process
of the hidden persistent component of the demand shock.

That scalar
innovation process is the additional state variable contributed by the
problem of extracting a signal from equilibrium prices that decision
makers face in Townsend’s model.






.. rubric:: Footnotes

.. [#footnote0] :cite:`Pearlman_Sargent2005` verified this assertion using a different tactic, namely, by constructing
      analytic formulas an equilibrium under the incomplete
      information structure and confirming that they match the pooling equilibrium formulas derived here.


.. [#footnote1]  See :cite:`ahms` for an account of invariant subspace methods.

.. [#footnote2]  See :cite:`ahms` for a discussion
      of the information assumptions needed to create a situation
      in which higher order beliefs appear in equilibrium decision rules.  The way
      to read our findings in light of :cite:`ams` is that Townsend's
      section 8 model  has too few sources of random shocks relative
      to sources of signals to permit higher order beliefs to
      play a role.

.. [#footnote3]  See :cite:`Sargent1987`, especially
      chapters IX and XIV, for the principles  that guide solving some roots backwards and others forwards.

.. [#footnote4]  As noted :cite:`Sargent1987`, this difference equation is the Euler equation for
      the planning problem   of maximizing the discounted sum of consumer plus
      producer surplus.

.. [#footnote5] :cite:`Pearlman_Sargent2005` verify the same claim by applying   machinery of  :cite:`PCL`.
