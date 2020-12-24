System 1
========

.. math::

    \begin{aligned}
        \left[\begin{array}{c}
        e_{t+1}\\
        k_{t+1}^{i}\\
        \tilde{\theta}_{t+1}\\
        P_{t+1}\\
        \theta_{t+1}\\
        v_{t+1}
        \end{array}\right] & = & \underbrace{\left[\begin{array}{cccccc}
        0 & 0 & 0 & 0 & 0 & 0\\
        \frac{1}{\lambda-\rho}\frac{\rho p}{p+\sigma_{e}^{2}} & \tilde{\lambda} & \frac{1}{\lambda-\rho}\frac{-\rho\sigma_{e}^{2}}{p+\sigma_{e}^{2}} & 0 & \frac{\rho}{\lambda-\rho} & 0\\
        -\frac{\rho p}{p+\sigma_{e}^{2}} & 0 & \frac{\rho\sigma_{e}^{2}}{p+\sigma_{e}^{2}} & 0 & 0 & 1\\
        \frac{b}{\lambda-\rho}\frac{\rho p}{p+\sigma_{e}^{2}} & b\tilde{\lambda} & b\frac{1}{\lambda-\rho}\frac{-\rho\sigma_{e}^{2}}{p+\sigma_{e}^{2}} & 0 & b\frac{\rho}{\lambda-\rho}+\rho & 1\\
        0 & 0 & 0 & 0 & \rho & 1\\
        0 & 0 & 0 & 0 & 0 & 0
        \end{array}\right]}_{A}\left[\begin{array}{c}
        e_{t}\\
        k_{t}^{i}\\
        \tilde{\theta}_{t}\\
        P_{t}\\
        \theta_{t}\\
        v_{t}
        \end{array}\right]+\underbrace{\left[\begin{array}{cc}
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
        G & = & \left[\begin{array}{cccccc}
        0 & 0 & 0 & 1 & 0 & 0\\
        1 & 0 & 0 & 0 & 1 & 0\\
        1 & 0 & 0 & 0 & 0 & 0
        \end{array}\right]\\
        H & = & \left[\begin{array}{c}
        0\\
        0\\
        0
        \end{array}\right]\\
        \left[\begin{array}{c}
        z_{1,t+1}\\
        z_{2,t+1}
        \end{array}\right] & \sim & \mathcal{N}\left(0,I\right)
    \end{aligned}

.. code-block:: ipython

    import numpy as np
    import quantecon as qe
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    import plotly.express as px
    import plotly.offline as pyo
    from statsmodels.regression.linear_model import OLS
    
    
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

.. code-block:: python3

    # Verify that p = σ_v ^ 2 + p * ρ ^ 2 - (ρ * p) ^ 2 / (p + σ_e ** 2)
    tol = 1e-12
    np.abs(p - (σ_v ** 2 + p * ρ ** 2 - (ρ * p) ** 2 / (p + σ_e ** 2))) < tol

.. code-block:: python3

    term_0 = -ρ * σ_e ** 2 / (p + σ_e ** 2)
    term_1 = ρ * p / (p + σ_e ** 2)
    
    A_lss = np.array([[0., 0., 0., 0., 0., 0.],
                     [term_1 / (λ - ρ), λ_tilde, term_0 / (λ - ρ), 0., ρ / (λ - ρ), 0.],
                     [-term_1, 0., -term_0, 0., 0., 1.],
                     [b * term_1 / (λ - ρ) , b * λ_tilde, b * term_0 / (λ - ρ), 0., b * ρ / (λ - ρ) + ρ, 1.],
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
    x, y = lss.simulate(ts_length)

.. code-block:: python3

    # Verify that two ways of computing P_t match
    np.max(np.abs(np.array([[1., b, 0., 0., 1., 0.]]) @ x - x[3])) < 1e-12

.. code-block:: python3

    # Plot sample time path
    t = 300
    
    subplot_titles = [r'$e_t$',
                      r'$k^{i}_t$',
                      r'$\tilde{\theta}_t$',
                      r'$P_t$',
                      r'$\theta_t$',
                     r'$v_t$']
    
    fig = make_subplots(rows=x.shape[0], cols=1, subplot_titles=subplot_titles)
    
    for idx in range(x.shape[0]):
        fig.add_trace(go.Scatter(y=x[idx, :t],
                             legendgroup='trend'),
                 row=idx+1,
                 col=1)
        
    fig.update_layout(height=1200)    
        
    fig.show()

.. code-block:: python3

    # Plot sample time path
    t = 300
    
    subplot_titles = [r'$P_t$',
                      r'$\theta_t + e_t$',
                      r'$e_t$']
    
    fig = make_subplots(rows=y.shape[0], cols=1, subplot_titles=subplot_titles)
    
    for idx in range(y.shape[0]):
        fig.add_trace(go.Scatter(y=y[idx, :t],
                             legendgroup='trend'),
                 row=idx+1,
                 col=1)
        
    fig.update_layout(height=800)    
        
    fig.show()

.. code-block:: python3

    fig = px.histogram(x[2], title=r'$\mathrm{Histogram: }\: \tilde{\theta}_{t}$')
    fig.update_layout(height=500)   

.. code-block:: python3

    # Compute the mean of \tilde{\theta}
    x[2].mean()

.. code-block:: python3

    xcoef, ycoef = lss.impulse_response(j=20)
    data = np.array([xcoef])[0, :, 1, :]
    
    fig = go.Figure(data=go.Scatter(y=data[:, 0], name=r'$e_{t+1}$'))
    fig.add_trace(go.Scatter(y=data[:, 1], name=r'$v_{t+1}$'))
    fig.update_layout(title=r'Impulse Response Function',
                       xaxis_title='Time',
                       yaxis_title=r'$k^{i}_{t}$')
    fig1=fig
    fig1.show()

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

System 2
========

.. math::

    \begin{aligned}
        \left[\begin{array}{c}
        e_{1,t+1}\\
        e_{2,t+1}\\
        k_{t+1}^{i}\\
        \tilde{\theta}_{t+1}\\
        P_{t+1}^{1}\\
        P_{t+1}^{2}\\
        \theta_{t+1}\\
        v_{t+1}
        \end{array}\right] & = & \underbrace{\left[\begin{array}{cccccccc}
        0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
        0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
        \frac{1}{\lambda-\rho}\frac{\rho p}{2p+\sigma_{e}^{2}} & \frac{1}{\lambda-\rho}\frac{\rho p}{2p+\sigma_{e}^{2}} & \tilde{\lambda} & \frac{1}{\lambda-\rho}\frac{-\rho\sigma_{e}^{2}}{2p+\sigma_{e}^{2}} & 0 & 0 & \frac{\rho}{\lambda-\rho} & 0\\
        -\frac{\rho p}{2p+\sigma_{e}^{2}} & -\frac{\rho p}{2p+\sigma_{e}^{2}} & 0 & \frac{\rho\sigma_{e}^{2}}{2p+\sigma_{e}^{2}} & 0 & 0 & 0 & 1\\
        \frac{b}{\lambda-\rho}\frac{\rho p}{2p+\sigma_{e}^{2}} & \frac{b}{\lambda-\rho}\frac{\rho p}{2p+\sigma_{e}^{2}} & b\tilde{\lambda} & b\frac{1}{\lambda-\rho}\frac{-\rho\sigma_{e}^{2}}{2p+\sigma_{e}^{2}} & 0 & 0 & b\frac{\rho}{\lambda-\rho}+\rho & 1\\
        \frac{b}{\lambda-\rho}\frac{\rho p}{2p+\sigma_{e}^{2}} & \frac{b}{\lambda-\rho}\frac{\rho p}{2p+\sigma_{e}^{2}} & b\tilde{\lambda} & b\frac{1}{\lambda-\rho}\frac{-\rho\sigma_{e}^{2}}{2p+\sigma_{e}^{2}} & 0 & 0 & b\frac{\rho}{\lambda-\rho}+\rho & 1\\
        0 & 0 & 0 & 0 & 0 & 0 & \rho & 1\\
        0 & 0 & 0 & 0 & 0 & 0 & 0 & 0
        \end{array}\right]}_{A}\left[\begin{array}{c}
        e_{1,t}\\
        e_{2,t}\\
        k_{t}^{i}\\
        \tilde{\theta}_{t}\\
        P_{t}^{1}\\
        P_{t}^{2}\\
        \theta_{t}\\
        v_{t}
        \end{array}\right]+\underbrace{\left[\begin{array}{ccc}
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
        G & = & \left[\begin{array}{cccccccc}
        0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\
        0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\\
        1 & 0 & 0 & 0 & 0 & 0 & 1 & 0\\
        0 & 1 & 0 & 0 & 0 & 0 & 1 & 0\\
        1 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
        0 & 1 & 0 & 0 & 0 & 0 & 0 & 0
        \end{array}\right]\\
        H & = & \left[\begin{array}{c}
        0\\
        0\\
        0\\
        0\\
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
    
    A_lss = np.array([[0., 0., 0., 0., 0., 0., 0., 0.],
                     [0., 0., 0., 0., 0., 0., 0., 0.],
                     [term_1 / (λ - ρ), term_1 / (λ - ρ), λ_tilde, term_0 / (λ - ρ), 0., 0., ρ / (λ - ρ), 0.],
                     [-term_1, -term_1, 0., -term_0, 0., 0., 0., 1.],
                     [b * term_1 / (λ - ρ), b * term_1 / (λ - ρ), b * λ_tilde, b * term_0 / (λ - ρ), 0., 0., b * ρ / (λ - ρ) + ρ, 1.],
                     [b * term_1 / (λ - ρ), b * term_1 / (λ - ρ), b * λ_tilde, b * term_0 / (λ - ρ), 0., 0., b * ρ / (λ - ρ) + ρ, 1.],
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
    x, y = lss.simulate(ts_length)

.. code-block:: python3

    # Plot sample time path
    t = 300
    
    subplot_titles = [r'$e_{1,t}$',
                      r'$e_{2,t}$',
                      r'$k^{i}_t$',
                      r'$\tilde{\theta}_t$',
                      r'$P^{1}_t$',
                      r'$P^{2}_t$',
                      r'$\theta_t$',
                     r'$v_t$']
    
    fig = make_subplots(rows=x.shape[0], cols=1, subplot_titles=subplot_titles)
    
    for idx in range(x.shape[0]):
        fig.add_trace(go.Scatter(y=x[idx, :t],
                             legendgroup='trend'),
                 row=idx+1,
                 col=1)
        
    fig.update_layout(height=1400)    
        
    fig.show()

.. code-block:: python3

    # Plot sample time path
    t = 300
    
    subplot_titles = [r'$P^{1}_t$',
                      r'$P^{2}_t$',
                      r'$\theta_t + e_{1,t}$',
                      r'$\theta_t + e_{2,t}$',
                      r'$e_{1,t}$',
                      r'$e_{2,t}$']
    
    fig = make_subplots(rows=y.shape[0], cols=1, subplot_titles=subplot_titles)
    
    for idx in range(y.shape[0]):
        fig.add_trace(go.Scatter(y=y[idx, :t],
                             legendgroup='trend'),
                 row=idx+1,
                 col=1)
        
    fig.update_layout(height=1200)    
        
    fig.show()

.. code-block:: python3

    xcoef, ycoef = lss.impulse_response(j=20)

.. code-block:: python3

    data = np.array([xcoef])[0, :, 2, :]
    
    fig = go.Figure(data=go.Scatter(y=data[:, 0], name=r'$e_{1,t+1}$'))
    fig.add_trace(go.Scatter(y=data[:, 1], name=r'$e_{2,t+1}$'))
    fig.add_trace(go.Scatter(y=data[:, 2], name=r'$v_{t+1}$'))
    fig.update_layout(title=r'Impulse Response Function',
                       xaxis_title='Time',
                       yaxis_title=r'$k^{i}_{t}$')
    fig2=fig
    fig2.show()


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

.. code-block:: python3

    # θ_t + e^{2}_t on k^{i}_t, P^{1}_t, P^{2}_t, \\tilde{\\theta_t}

.. code-block:: python3

    # Verify that θ_t + e^{2}_t can be recovered
    model = OLS(y[1], x[2:6].T)
    reg_res = model.fit()
    np.abs(reg_res.rsquared - 1.) < 1e-6 

.. code-block:: python3

    reg_res.rsquared

.. code-block:: python3

    fig1.show()

.. code-block:: python3

    fig2.show()

