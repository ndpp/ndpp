.. _methods_chi:

========================
Chi Integration Overview
========================

*NOTE THAT THIS FEATURE HAS NOT YET BEEN TESTED*
Fissionable nuclides can contain one ore more fission reactions in their ACE
data.  NDPP must parse through each of these reactions (and the associated
delayed neutron precursor data), calculate the outgoing
energy spectra, :math:`\chi` for each energy group at each incoming energy point
in the ACE tables (:math:`\chi_g(E')`), and then combine these in to a single
value of :math:`\chi_g(E')` (on a single energy grid). To generate
:math:`\chi_g(E')`, NDPP is solving the following equations for the prompt,
delayed, and total values of :math:`\chi_g(E')`:

.. math::

   \chi_{prompt,g}(E') =\ \sum\limits_{MT}\:\sum\limits_{d}
       \frac{\sigma_d(E')}{\sigma_f(E')}\int\limits_{E_g}^{E_{g-1}}
       \chi_{MT,d}\left(E,E'\right)dE'

.. math::

   \chi_{delayed,g}(E') =\ \sum\limits_{c}\:\ Y_c(E')
       \int\limits_{E_g}^{E_{g-1}}\chi_c\left(E,E'\right)dE'

.. math::

   \chi_{total,g}(E') =\ \left(1-\beta(E')\right) \chi_{prompt,g}(E') +
       \beta(E') \chi_{delayed,g}(E')

In the above equations, :math:`E'` is the incoming neutron energy, :math:`MT`
is the reaction channel, :math:`d` is the energy distributions within that
channel, :math:`E_g` and :math:`E_{g-1}` are the lower and upper energy group
boundaries for group :math:`g` respectively, :math:`\beta(E')` is the total
delayed neutron emission fraction, :math:`\sigma_{d,MT}(E')` and
:math:`\sigma_{f,tot}(E')` are  the microscopic cross-sections of this reaction
channel and distribution occuring and the microscopic cross-section for all
fission reactions, and :math:`Y_c` is the yield of precursor group :math:`c`.

The details of this process are discussed in the following sections.

------------------------------------
Fission Reactions and Initialization
------------------------------------

To determine if a nuclide is fissionable, NDPP checks for presence of the
:math:`MT=18` reaction channel in the ACE data. This channel is the total
fission reaction, and must be present for there to be any possibility of
generating values of :math:`\chi_g(E')` for the nuclide. If the nuclide is
fissionable, then the code determines the number of fission reaction channels,
outgoing energy distributions, delayed neutron precursor groups, and the number
of energy grid points for each.

Then, NDPP progresses through each fission channel (and energy distributions
within that channel) and delayed neutron precursor group to obtain values for
each of the terms in the :math:`\chi_g(E')` equation above.  The values of
:math:`\beta`, :math:`\sigma_d`, :math:`\sigma_f`, and :math:`Y_c` are
calculated by the same means discussed in the OpenMC manual. The calculation of
:math:`\int\limits_{E_g}^{E_{g-1}}\chi\left(E,E'\right)dE'` is described in the
next section.

-------------------------------------
Outgoing Neutron Energy Distributions
-------------------------------------

Monte Carlo codes must sample single values from probability distribution
functions; NDPP, however, must integrate that probability distribution function
between upper and lower boundaries.  Therefore, the ACE Laws are utilized in
different ways than described in the OpenMC manual and will be discussed herein.

All fission reactions are obtained from either ENDF File 5 (Energy Distribution
of Secondary Particles) or File 6 (Product Energy-Angle Distributions).  The
following subsections discuss each of the ACE Laws used to describe fission
neutron energy distributions in ENDF/B-7 and how these laws are treated by NDPP.

After the calculation of
:math:`\int\limits_{E_g}^{E_{g-1}}\chi\left(E,E'\right)dE'` for each energy
group, the values are normalized to 1.0 to account for any inaccuracies
introduced by the interpolation schemes.

ACE Law 4 - Continuous Tabular Distribution
+++++++++++++++++++++++++++++++++++++++++++

This representation is essentially a two-dimensional table which provides
points of both the probability distribution function
:math:`\chi\left(E,E'\right)`, a cumulative distribution function (CDF)
:math:`\int\limits_{0}^E\chi\left(E,E'\right)dE'`, and rules for interpolating
between each of the E and E' data sets.  To determine
:math:`\int\limits_{E_g}^{E_{g-1}}\chi\left(E,E'\right)dE'`, NDPP must first
find the location of the data corresponding to :math:`E'` data, and then
interpolate on the CDF to find teh CDF at :math:`E_{g-1}` and :math:`E_g`.  The
value of :math:`\int\limits_{E_g}^{E_{g-1}}\chi\left(E,E'\right)dE'` is simply
the difference between these two values.  The interpolation rules are followed
as described in the OpenMC Methods Manual.

ACE Law 7 - Maxwell Fission Spectrum
++++++++++++++++++++++++++++++++++++

One representation of the secondary energies for neutrons from fission is the
so-called Maxwell spectrum. A probability distribution for the Maxwell spectrum
can be written in the form

.. math::
    \chi(E,E') dE' = c E'^{1/2} e^{-E'/T(E)} dE'

where :math:`E` is the incoming energy of the neutron and :math:`T` is the
so-called nuclear temperature, which is a function of the incoming energy of the
neutron. The ACE format contains a list of nuclear temperatures versus incoming
energies. The nuclear temperature is interpolated between neighboring incoming
energies using a specified interpolation law. Once the temperature :math:`T` is
determined, we then can analytically determine the value of
:math:`\int\limits_{E_g}^{E_{g-1}}\chi\left(E',E\right)dE` with the following
relation:

.. math::
    \int\limits_{E_g}^{E_{g-1}}\chi\left(E,E'\right)dE' =\
        c \left(\frac{1}{2}\sqrt{\pi}\left(T(E)\right)^{\frac{3}{2}}
        erf\left(\frac{E'}{T(E)}\right)-T(E)\sqrt{E'}\exp{-\frac{E'}{T(E)}}\right)

This integral is forced to 0 for values of E' greater than the restriction
energy, :math:`U(E)`.

ACE Law 9 - Evaporation Spectrum
++++++++++++++++++++++++++++++++

Evaporation spectra are primarily used in compound nucleus processes where a
secondary particle can "evaporate" from the compound nucleus if it has
sufficient energy. The probability distribution for an evaporation spectrum can
be written in the form

.. math::
    \chi(E,E') dE' = c E' e^{-E'/T(E)} dE'

where :math:`E` is the incoming energy of the neutron and :math:`T` is the
nuclear temperature, which is a function of the incoming energy of the
neutron. The ACE format contains a list of nuclear temperatures versus incoming
energies. The nuclear temperature is interpolated between neighboring incoming
energies using a specified interpolation law. Once the temperature :math:`T` is
determined, we then analytically determine the value of
:math:`\int\limits_{E_g}^{E_{g-1}}\chi\left(E,E'\right)dE'` with the following
relation:

.. math::
    \int\limits_{E_g}^{E_{g-1}}\chi\left(E,E'\right)dE' =\
        -T(E) c \exp{-\frac{E'}{T(E)}}\left(T(E)+E'\right)

This integral is forced to 0 for values of E' greater than the restriction
energy, :math:`U(E)`.

ACE Law 11 - Energy-Dependent Watt Spectrum
+++++++++++++++++++++++++++++++++++++++++++

The probability distribution for a Watt fission spectrum can be written in the
form

.. math::
    \chi(E,E') dE' = c e^{-E'/a(E)} \sinh \sqrt{b(E) \, E'} dE'

where :math:`a` and :math:`b` are parameters for the distribution and are given
as tabulated functions of the incoming energy of the neutron. These two
parameters are interpolated on the incoming energy grid using a specified
interpolation law. Once the parameters have been determined, we then
analytically determine the value of
:math:`\int\limits_{E_g}^{E_{g-1}}\chi\left(E,E'\right)dE'` with the following
relation:

.. .. math::
    \int\limits_{E_g}^{E_{g-1}}\chi\left(E,E'\right)dE' =\

This integral is forced to 0 for values of E' greater than the restriction
energy, :math:`U(E)`.

ACE Law 61 - Correlated Energy and Angle Distribution
+++++++++++++++++++++++++++++++++++++++++++++++++++++

This law is very similar to ACE Law 4, except there is another dimension in the
table to represent the angular probability distribution function.  Since the
:math:`\chi` portion of NDPP is not concerned with the outgoing angle, and
therefore this extra dimension can be ignored.  Therefore the methods used to
calculate :math:`\int\limits_{E_g}^{E_{g-1}}\chi\left(E,E'\right)dE'`, are the
same as is discussed in the Law 4 section.

----------------------------------------------------
Creation of Union Energy Grids for :math:`\chi_g(E)`
----------------------------------------------------

At this stage, NDPP has a tabular representation of
:math:`\int\limits_{E_g}^{E_{g-1}}\chi\left(E,E'\right)dE'` for each incoming
energy, :math:`E`, and outgoing energy group, :math:`g` for every fission
reaction channel and energy distribution as well as for each of the delayed
neutron precursor groups.  Each of these tables has values on a completely
different set of incoming energies (since the ACE data are on separate energy
grids as well) and must be combined on to the same energy grid for the prompt,
delayed, and total values of :math:`\chi_g(E)`.  This unionized energy grid
is made by using all of the energy points in the relevant :math:`\chi_g(E)`
distributions and linearly interpolating between values for points without a
data set on the grid. Due to the additional interpolation step, these values
are also re-normalized to 1.0.  A unioninzed grid exists for each of the prompt,
delayed, and total values of :math:`\chi_g(E)`.

------------------------------
Thinning of Union Energy Grids
------------------------------

Since the unionized grids must be searched by the Monte Carlo code during
runtime, it is desirable to have the size of the grid be as small as possible.
To this end, NDPP provides the user with an option to `thin` the energy grid
such that :math:`E` points which provide an increase accuracy of less than the
user-specified tolerance when linear interpolation with neighboring points is
used instead of the explicit value are discared from the data.



