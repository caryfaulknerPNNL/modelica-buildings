within Buildings;
package AirCleaning

  model PAC "In-room portable air cleaner"

    replaceable package Medium =
      Modelica.Media.Interfaces.PartialMedium "Medium in the component";

    parameter Real eff[Medium.nC](min=0, max=1) = 0.9997
      "Virus removal efficiency";

    parameter Integer nPACs(min=0) = 1
      "Number of PACs";

    parameter Real flowPAC(min=0) = 0.094
      "PAC flow rate";

    parameter Real kpow(min=0) = 50
      "Rated power";

    Modelica.Blocks.Interfaces.RealInput C[Medium.nC] "Zone concentration"
      annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
    Modelica.Blocks.Interfaces.RealOutput yC_flow[Medium.nC] "Concentration outflow"
      annotation (Placement(transformation(extent={{100,30},{120,50}})));
    Modelica.Blocks.Interfaces.BooleanInput uPACEna "True when PAC is on"
      annotation (Placement(transformation(extent={{-140,-40},{-100,0}})));
    Modelica.Blocks.Interfaces.RealOutput yP_PAC "Power output"
      annotation (Placement(transformation(extent={{100,-10},{120,10}})));
    Modelica.Blocks.Interfaces.RealOutput yE_PAC "Energy output"
      annotation (Placement(transformation(extent={{100,-50},{120,-30}})));
    Modelica.Blocks.Math.BooleanToReal booleanToReal
      annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
    Modelica.Blocks.Continuous.Integrator E_PAC(use_reset=false)
      annotation (Placement(transformation(extent={{60,-50},{80,-30}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a diss "Heat dissipation"
      annotation (Placement(transformation(extent={{94,-90},{114,-70}})));

  protected
    Modelica.Blocks.Math.Gain pow(final k=nPACs*kpow) "power of PAC"
      annotation (Placement(transformation(extent={{-20,-30},{0,-10}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prePow(final alpha=0)
      "Prescribed power (=heat and flow work) flow for dynamic model"
      annotation (Placement(transformation(extent={{10,-70},{30,-50}})));

  equation
    connect(uPACEna, booleanToReal.u)
      annotation (Line(points={{-120,-20},{-82,-20}}, color={255,0,255}));
    connect(booleanToReal.y, pow.u)
      annotation (Line(points={{-59,-20},{-22,-20}}, color={0,0,127}));
    connect(pow.y,yP_PAC)  annotation (Line(points={{1,-20},{96,-20},{96,0},{110,
            0}}, color={0,0,127}));
    connect(pow.y, prePow.Q_flow) annotation (Line(points={{1,-20},{8,-20},{8,-46},
            {2,-46},{2,-60},{10,-60}}, color={0,0,127}));
    connect(prePow.port, diss) annotation (Line(points={{30,-60},{90,-60},{90,-80},
            {104,-80}}, color={191,0,0}));
    connect(pow.y,E_PAC. u) annotation (Line(points={{1,-20},{8,-20},{8,-40},{58,
            -40}}, color={0,0,127}));
    connect(E_PAC.y,yE_PAC)
      annotation (Line(points={{81,-40},{110,-40}}, color={0,0,127}));

    for i in 1:Medium.nC loop
      yC_flow[i]  = booleanToReal.y*(-flowPAC)*eff[i]*nPACs*C[i];
    end for;
      annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"),
                  Dialog(group = "Integer"),
                Icon(graphics={
          Ellipse(
            extent={{-30,68},{32,50}},
            lineColor={28,108,200},
            fillColor={28,108,200},
            fillPattern=FillPattern.None),
          Ellipse(
            extent={{-30,-42},{32,-60}},
            lineColor={28,108,200},
            fillColor={28,108,200},
            fillPattern=FillPattern.None),
          Line(points={{-30,60},{-30,-50}}, color={28,108,200}),
          Line(points={{32,60},{32,-50}}, color={28,108,200})}), Documentation(
          info="<html>
<p>This model is designed to simulate the removal of airborne contaminants using a portable air cleaner(PAC).The PAC model calculates the mass flow rate of the contaminants of the outflow from the PAC(<span style=\"font-family: Courier New;\">yC_flow[]</span>) based on the equation described below when the PAC is proven on(<span style=\"font-family: Courier New;\">uPACEna</span>). The PACs are modeled to consume energy using a constant power rating (<code>kpow</code>) and heat is dissipated into the zone based on the consumed power(<span style=\"font-family: Courier New;\">yP_PAC</span>).<span style=\"font-family: Courier New;\">yP_PAC</span>, and total energy output (<span style=\"font-family: Courier New;\">yE_PAC</span>), are calculated using the number of PACs in the zone(<span style=\"font-family: Courier New;\">nPACs</span>) and (<span style=\"font-family: Courier New;\">kpow</span>). </p>
<h4>Main Equations </h4>
<p align=\"center\"><span style=\"font-family: Courier New;\">Ċ</span><sub>PAC</sub> = eff<sub>PAC</sub>* <span style=\"font-family: Courier New;\">V̇</span><sub>PAC</sub>*c<sub>zone</sub> </p>
<p>where <span style=\"font-family: Courier New;\">Ċ</span><sub>PAC</sub> is the rate of pathogen removal by the PAC, eff<sub>PAC</sub> the pathogen removal efficency, <span style=\"font-family: Courier New;\">V̇</span><sub>PAC</sub> is the volumetric airflow rate of the PAC, and c<sub>zone</sub> is the pathogen concentration in the zone where the PAC is located.</p>
<h4>Assumptions</h4>
The model assumes well-mixed zones with uniform concentrations.

</html>"));
  end PAC;

  model InDuctGUV "In Duct GUV"
    extends Buildings.Fluid.BaseClasses.PartialInDuctGUV(
      final m_flow_turbulent = if computeFlowResistance then deltaM * m_flow_nominal_pos else 0, vol(
          nPorts=2));

    parameter Real deltaM(min=1E-6) = 0.3
      "Fraction of nominal mass flow rate where transition to turbulent occurs"
         annotation(Evaluate=true,
                    Dialog(group = "Transition to laminar",
                           enable = not linearized));

    parameter Real kpow(min=0) = 120
      "Rated power";

    parameter Boolean addPowerToMedium=true
      "Set to false to avoid any power (=heat and flow work) being added to medium (may give simpler equations)";

    final parameter Real k = if computeFlowResistance then
          m_flow_nominal_pos / sqrt(dp_nominal_pos) else 0
      "Flow coefficient, k=m_flow/sqrt(dp), with unit=(kg.m)^(1/2)";
    Modelica.Blocks.Interfaces.BooleanInput u "on/off"
      annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));
    Modelica.Blocks.Math.BooleanToReal booleanToReal
      annotation (Placement(transformation(extent={{-60,-90},{-40,-70}})));
    Modelica.Blocks.Continuous.Integrator E_GUV(use_reset=false)
      annotation (Placement(transformation(extent={{40,-80},{60,-60}})));
    Modelica.Blocks.Interfaces.RealOutput yP_GUV "Power output"
      annotation (Placement(transformation(extent={{100,-50},{120,-30}})));
    Modelica.Blocks.Interfaces.RealOutput yE_GUV "Energy output"
      annotation (Placement(transformation(extent={{100,-90},{120,-70}})));
    HVACFilter                        filt(
      allowFlowReversal=true,
      eff=0.87,
      m_flow_nominal=m_flow_nominal,
      dp_nominal=0,
      redeclare package Medium = Medium,
      use_eff=true)
      annotation (Placement(transformation(extent={{50,-10},{70,10}})));
    Modelica.Blocks.Sources.RealExpression mFlowRat(y=booleanToReal.y*(1 - exp(-2.18
          *m_flow_nominal/2/m_flow))) "flow ratio kinda"
      annotation (Placement(transformation(extent={{-72,22},{-52,42}})));
  protected
    final parameter Boolean computeFlowResistance=(dp_nominal_pos > Modelica.Constants.eps)
      "Flag to enable/disable computation of flow resistance"
     annotation(Evaluate=true);
    final parameter Real coeff=
      if linearized and computeFlowResistance
      then if from_dp then k^2/m_flow_nominal_pos else m_flow_nominal_pos/k^2
      else 0
      "Precomputed coefficient to avoid division by parameter";
  protected
    Modelica.Blocks.Math.Gain pow(final k=kpow)
                                               "power of GUV"
      annotation (Placement(transformation(extent={{-20,-60},{0,-40}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prePow(final alpha=0) if
      addPowerToMedium
      "Prescribed power (=heat and flow work) flow for dynamic model"
      annotation (Placement(transformation(extent={{-2,-90},{18,-70}})));
  protected
    Modelica.Blocks.Math.Gain keff(final k=eff)
                                               "power of GUV"
      annotation (Placement(transformation(extent={{-18,28},{2,48}})));
  initial equation
   if computeFlowResistance then
     assert(m_flow_turbulent > 0, "m_flow_turbulent must be bigger than zero.");
   end if;

   assert(m_flow_nominal_pos > 0, "m_flow_nominal_pos must be non-zero. Check parameters.");
  equation
    // Pressure drop calculation

    /*if computeFlowResistance then
    if linearized then
      if from_dp then
        m_flow = dp*coeff;
      else
        dp = m_flow*coeff;
      end if;
    else
      if homotopyInitialization then
        if from_dp then
          m_flow=homotopy(
            actual=Buildings.Fluid.BaseClasses.FlowModels.basicFlowFunction_dp(
              dp=dp,
              k=k,
              m_flow_turbulent=m_flow_turbulent),
            simplified=m_flow_nominal_pos*dp/dp_nominal_pos);
        else
          dp=homotopy(
            actual=Buildings.Fluid.BaseClasses.FlowModels.basicFlowFunction_m_flow(
              m_flow=m_flow,
              k=k,
              m_flow_turbulent=m_flow_turbulent),
            simplified=dp_nominal_pos*m_flow/m_flow_nominal_pos);
         end if;  // from_dp
      else // do not use homotopy
        if from_dp then
          m_flow=Buildings.Fluid.BaseClasses.FlowModels.basicFlowFunction_dp(
            dp=dp,
            k=k,
            m_flow_turbulent=m_flow_turbulent);
        else
          dp=Buildings.Fluid.BaseClasses.FlowModels.basicFlowFunction_m_flow(
            m_flow=m_flow,
            k=k,
            m_flow_turbulent=m_flow_turbulent);
        end if;  // from_dp
      end if; // homotopyInitialization
    end if; // linearized
  else // do not compute flow resistance
    dp = 0;
    end if;  // computeFlowResistance */

    connect(u, booleanToReal.u)
      annotation (Line(points={{-120,-80},{-62,-80}}, color={255,0,255}));
    connect(booleanToReal.y, pow.u) annotation (Line(points={{-39,-80},{-28,-80},{
            -28,-50},{-22,-50}}, color={0,0,127}));
    connect(pow.y, E_GUV.u) annotation (Line(points={{1,-50},{32,-50},{32,-70},{38,
            -70}}, color={0,0,127}));
    connect(pow.y, yP_GUV) annotation (Line(points={{1,-50},{32,-50},{32,-36},{94,
            -36},{94,-40},{110,-40}}, color={0,0,127}));
    connect(E_GUV.y, yE_GUV) annotation (Line(points={{61,-70},{94,-70},{94,-80},{
            110,-80}}, color={0,0,127}));
    connect(pow.y, prePow.Q_flow) annotation (Line(points={{1,-50},{6,-50},{6,-66},
            {-10,-66},{-10,-80},{-2,-80}}, color={0,0,127}));
    connect(port_a, vol.ports[1]) annotation (Line(points={{-100,0},{-16,0},{-16,-14},
            {-2,-14},{-2,-10}}, color={0,127,255}));
    connect(prePow.port, vol.heatPort) annotation (Line(points={{18,-80},{22,-80},
            {22,16},{-14,16},{-14,0},{-10,0}}, color={191,0,0}));
    connect(vol.ports[2], filt.port_a) annotation (Line(points={{2,-10},{2,-12},{-16,
            -12},{-16,18},{44,18},{44,0},{50,0}}, color={0,127,255}));
    connect(filt.port_b, port_b)
      annotation (Line(points={{70,0},{100,0}}, color={0,127,255}));
    connect(keff.y, filt.eff1) annotation (Line(points={{3,38},{22,38},{22,16},{
            42,16},{42,-6},{48,-6}}, color={0,0,127}));
    connect(mFlowRat.y, keff.u) annotation (Line(points={{-51,32},{-28,32},{-28,
            38},{-20,38}}, color={0,0,127}));
    annotation (defaultComponentName="res",
  Documentation(info="<html>
<p>
Model of a flow resistance with a fixed flow coefficient.
The mass flow rate is
</p>
<p align=\"center\" style=\"font-style:italic;\">
m&#775; = k
&radic;<span style=\"text-decoration:overline;\">&Delta;p</span>,
</p>
<p>
where
<i>k</i> is a constant and
<i>&Delta;p</i> is the pressure drop.
The constant <i>k</i> is equal to
<code>k=m_flow_nominal/sqrt(dp_nominal)</code>,
where <code>m_flow_nominal</code> and <code>dp_nominal</code>
are parameters.
</p>
<h4>Assumptions</h4>
<p>
In the region
<code>abs(m_flow) &lt; m_flow_turbulent</code>,
the square root is replaced by a differentiable function
with finite slope.
The value of <code>m_flow_turbulent</code> is
computed as
<code>m_flow_turbulent = deltaM * abs(m_flow_nominal)</code>,
where <code>deltaM=0.3</code> and
<code>m_flow_nominal</code> are parameters that can be set by the user.
</p>
<p>
The figure below shows the pressure drop for the parameters
<code>m_flow_nominal=5</code> kg/s,
<code>dp_nominal=10</code> Pa and
<code>deltaM=0.3</code>.
</p>
<p align=\"center\">
<img alt=\"image\" src=\"modelica://Buildings/Resources/Images/Fluid/FixedResistances/PressureDrop.png\"/>
</p>
<h4>Important parameters</h4>
<p>
The parameter <code>from_dp</code> is used to determine
whether the mass flow rate is computed as a function of the
pressure drop (if <code>from_dp=true</code>), or vice versa.
This setting can affect the size of the nonlinear system of equations.
</p>
<p>
If the parameter <code>linearized</code> is set to <code>true</code>,
then the pressure drop is computed as a linear function of the
mass flow rate.
</p>
<p>
Setting <code>allowFlowReversal=false</code> can lead to simpler
equations. However, this should only be set to <code>false</code>
if one can guarantee that the flow never reverses its direction.
This can be difficult to guarantee, as pressure imbalance after
the initialization, or due to medium expansion and contraction,
can lead to reverse flow.
</p>
<p>
If the parameter
<code>show_T</code> is set to <code>true</code>,
then the model will compute the
temperature at its ports. Note that this can lead to state events
when the mass flow rate approaches zero,
which can increase computing time.
</p>
<h4>Notes</h4>
<p>
For more detailed models that compute the actual flow friction,
models from the package
<a href=\"modelica://Modelica.Fluid\">
Modelica.Fluid</a>
can be used and combined with models from the
<code>Buildings</code> library.
</p>
<p>
For a model that uses the hydraulic parameter and flow velocity at nominal conditions
as a parameter, use
<a href=\"modelica://Buildings.Fluid.FixedResistances.HydraulicDiameter\">
Buildings.Fluid.FixedResistances.HydraulicDiameter</a>.
</p>
<h4>Implementation</h4>
<p>
The pressure drop is computed by calling a function in the package
<a href=\"modelica://Buildings.Fluid.BaseClasses.FlowModels\">
Buildings.Fluid.BaseClasses.FlowModels</a>,
This package contains regularized implementations of the equation
</p>
<p align=\"center\" style=\"font-style:italic;\">
  m&#775; = sign(&Delta;p) k  &radic;<span style=\"text-decoration:overline;\">&nbsp;&Delta;p &nbsp;</span>
</p>
<p>
and its inverse function.
</p>
<p>
To decouple the energy equation from the mass equations,
the pressure drop is a function of the mass flow rate,
and not the volume flow rate.
This leads to simpler equations.
</p>
</html>",   revisions="<html>
<ul>
<li>
September 21, 2018, by Michael Wetter:<br/>
Decrease value of <code>deltaM(min=...)</code> attribute.
See <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1026\">#1026</a>.
</li>
<li>
February 3, 2018, by Filip Jorissen:<br/>
Revised implementation of pressure drop equation
such that it depends on <code>from_dp</code>
when <code>linearized=true</code>.
See <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/884\">#884</a>.
</li>
<li>
December 1, 2016, by Michael Wetter:<br/>
Simplified model by removing the geometry dependent parameters into the new
model
<a href=\"modelica://Buildings.Fluid.FixedResistances.HydraulicDiameter\">
Buildings.Fluid.FixedResistances.HydraulicDiameter</a>.
</li>
<li>
November 23, 2016, by Filip Jorissen:<br/>
Removed <code>dp_nominal</code> and
<code>m_flow_nominal</code> labels from icon.
</li>
<li>
October 14, 2016, by Michael Wetter:<br/>
Updated comment for parameter <code>use_dh</code>.
</li>
<li>
November 26, 2014, by Michael Wetter:<br/>
Added the required <code>annotation(Evaluate=true)</code> so
that the system of nonlinear equations in
<a href=\"modelica://Buildings.Fluid.FixedResistances.Validation.PressureDropsExplicit\">
Buildings.Fluid.FixedResistances.Validation.PressureDropsExplicit</a>
remains the same.
</li>
<li>
November 20, 2014, by Michael Wetter:<br/>
Rewrote the warning message using an <code>assert</code> with
<code>AssertionLevel.warning</code>
as this is the proper way to write warnings in Modelica.
</li>
<li>
August 5, 2014, by Michael Wetter:<br/>
Corrected error in documentation of computation of <code>k</code>.
</li>
<li>
May 29, 2014, by Michael Wetter:<br/>
Removed undesirable annotation <code>Evaluate=true</code>.
</li>
<li>
October 8, 2013, by Michael Wetter:<br/>
Removed parameter <code>show_V_flow</code>.
</li>
<li>
December 14, 2012 by Michael Wetter:<br/>
Renamed protected parameters for consistency with the naming conventions.
</li>
<li>
January 16, 2012 by Michael Wetter:<br/>
To simplify object inheritance tree, revised base classes
<code>Buildings.Fluid.BaseClasses.PartialResistance</code>,
<code>Buildings.Fluid.Actuators.BaseClasses.PartialTwoWayValve</code>,
<code>Buildings.Fluid.Actuators.BaseClasses.PartialDamperExponential</code>,
<code>Buildings.Fluid.Actuators.BaseClasses.PartialActuator</code>
and model
<code>Buildings.Fluid.FixedResistances.PressureDrop</code>.
</li>
<li>
May 30, 2008 by Michael Wetter:<br/>
Added parameters <code>use_dh</code> and <code>deltaM</code> for easier parameterization.
</li>
<li>
July 20, 2007 by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"),
      Icon(graphics={
          Rectangle(
            extent={{-74,94},{72,76}},
            lineColor={28,108,200},
            fillColor={85,170,255},
            fillPattern=FillPattern.Solid),
          Line(points={{-60,70},{-60,-14}}, color={28,108,200}),
          Line(points={{-40,70},{-40,-14}}, color={28,108,200}),
          Line(points={{0,70},{0,-14}}, color={28,108,200}),
          Line(points={{60,70},{60,-14}}, color={28,108,200}),
          Line(points={{40,70},{40,-14}}, color={28,108,200}),
          Line(points={{20,70},{20,-14}}, color={28,108,200}),
          Line(points={{-20,70},{-20,-14}}, color={28,108,200})}),
      experiment(StopTime=7200, __Dymola_Algorithm="Dassl"));
  end InDuctGUV;

  model GUV "In-room GUV"

    parameter Real frad(min=0, max=1) = 0.2
      "Fraction of irradiated space";

    parameter Real Eavg(min=0, max=1) = 50e-6
      "Effluence rate";

    parameter Real krad(min=0) = 2.93e3
      "Inactivation constant";

    parameter Real kpow(min=0) = 120
      "Rated power";

    parameter Real V(min=0) = 120
      "Zone volume";

    Modelica.Blocks.Interfaces.RealInput C "Zone concentration"
      annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
    Modelica.Blocks.Interfaces.RealOutput yC_flow "Concentration outflow"
      annotation (Placement(transformation(extent={{100,30},{120,50}})));
    Modelica.Blocks.Interfaces.BooleanInput u "on/off"
      annotation (Placement(transformation(extent={{-140,-40},{-100,0}})));
    Modelica.Blocks.Interfaces.RealOutput yP_GUV "Power output"
      annotation (Placement(transformation(extent={{100,-10},{120,10}})));
    Modelica.Blocks.Interfaces.RealOutput yE_GUV "Energy output"
      annotation (Placement(transformation(extent={{100,-50},{120,-30}})));
    Modelica.Blocks.Math.BooleanToReal booleanToReal
      annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
    Modelica.Blocks.Continuous.Integrator E_GUV(use_reset=false)
      annotation (Placement(transformation(extent={{60,-50},{80,-30}})));
    Modelica.Blocks.Math.MultiProduct multiProduct(nu=2)
      annotation (Placement(transformation(extent={{-80,34},{-68,46}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a diss "Heat dissipation"
      annotation (Placement(transformation(extent={{94,-90},{114,-70}})));

  protected
    Modelica.Blocks.Math.Gain pow(final k=kpow)
                                               "power of GUV"
      annotation (Placement(transformation(extent={{-20,-30},{0,-10}})));
    Modelica.Blocks.Math.Gain UV_calc(final k=-1.2*V*frad*Eavg*krad)
      "death of virus from UV"
      annotation (Placement(transformation(extent={{-16,30},{4,50}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prePow(final alpha=0)
      "Prescribed power (=heat and flow work) flow for dynamic model"
      annotation (Placement(transformation(extent={{10,-70},{30,-50}})));

  equation
    connect(u, booleanToReal.u)
      annotation (Line(points={{-120,-20},{-82,-20}}, color={255,0,255}));
    connect(C, multiProduct.u[1]) annotation (Line(points={{-120,40},{-100,40},{
            -100,42.1},{-80,42.1}}, color={0,0,127}));
    connect(booleanToReal.y, multiProduct.u[2]) annotation (Line(points={{-59,-20},
            {-54,-20},{-54,28},{-84,28},{-84,37.9},{-80,37.9}}, color={0,0,127}));
    connect(multiProduct.y, UV_calc.u)
      annotation (Line(points={{-66.98,40},{-18,40}}, color={0,0,127}));
    connect(booleanToReal.y, pow.u)
      annotation (Line(points={{-59,-20},{-22,-20}}, color={0,0,127}));
    connect(pow.y, yP_GUV) annotation (Line(points={{1,-20},{96,-20},{96,0},{110,
            0}}, color={0,0,127}));
    connect(pow.y, prePow.Q_flow) annotation (Line(points={{1,-20},{8,-20},{8,-46},
            {2,-46},{2,-60},{10,-60}}, color={0,0,127}));
    connect(prePow.port, diss) annotation (Line(points={{30,-60},{90,-60},{90,-80},
            {104,-80}}, color={191,0,0}));
    connect(pow.y, E_GUV.u) annotation (Line(points={{1,-20},{8,-20},{8,-40},{58,
            -40}}, color={0,0,127}));
    connect(E_GUV.y, yE_GUV)
      annotation (Line(points={{81,-40},{110,-40}}, color={0,0,127}));
    connect(UV_calc.y, yC_flow)
      annotation (Line(points={{5,40},{110,40}}, color={0,0,127}));
    annotation (Icon(graphics={
          Polygon(
            points={{-60,58},{-60,-64},{-40,-64},{-40,58},{-60,58}},
            lineColor={28,108,200},
            fillColor={28,108,200},
            fillPattern=FillPattern.Solid),
          Line(points={{-36,48},{42,48}}, color={28,108,200}),
          Line(points={{-36,26},{42,26}}, color={28,108,200}),
          Line(points={{-36,8},{42,8}},   color={28,108,200}),
          Line(points={{-36,-32},{42,-32}}, color={28,108,200}),
          Line(points={{-36,-52},{42,-52}}, color={28,108,200}),
          Line(points={{-36,-12},{42,-12}}, color={28,108,200})}),
        Documentation(info="<html>

<p>This model is designed to simulate the inactivation of airborne pathogens by either upper-or whole-room Germicidal Ultraviolet light devices(GUV) when the GUV is enabled <code>uGUVEna</code>. The GUV model calculates the mass flow rate of the active pathogens in the outflow from the GUV(<span style=\"font-family: Courier New;\">yC_flow</span>) based on the equation described below. The GUV device consumes energy based on a constant power rating,<code>kpow</code>, which is ultimately dissipated into the zone as heat(<code>diss</code>).
(<span style=\"font-family: Courier New;\">Kpow</span>) is used to calculate the power output(<span style=\"font-family: Courier New;\">yP_PAC</span>), and total energy output (<span style=\"font-family: Courier New;\">yE_PAC</span>). </p>
</p>

<h4> Main Equations </h4>
<p align=\"center\"><span style=\"font-family: Courier New;\">Ċ</span><sub>GUV</sub> = E<sub>avg</sub> <span style=\"font-family: Courier New;\">k</span><sub>rad</sub>f<sub>rad</sub>V<sub>zone</sub>c<sub>zone</sub> </p>
<p>where <span style=\"font-family: Courier New;\">Ċ</span><sub>PAC</sub> is the inactivation rate of airborne pathogens by the GUV, E<sub>avg</sub> is the average fluence rate of the GUV device, k</span><sub>rad</sub> is susceptibility of the pathogens to the GUV irrdiation,f<sub>rad</sub> is the fraction of irradiated volume in the zone, V<sub>zone</sub> is the total volume of the zone, and c<sub>zone</sub> is the pathogen concentration in the zone where the PAC is located.</p>

<h4> Assumptions </h4>

This assumes that the zone is well mixed and the parameters E<sub>avg</sub>,k<sub>rad</sub>, and f<sub>rad</sub> are constant for the given simulation when the GUV device is on.

</html>
"));
  end GUV;

  package BaseClasses
    model OutputVariableSelection
      "This model set the variables that will be output in .mat file"
      annotation (__Dymola_selections={Selection(name="SelectedVariable",match={
      MatchVariable(name="flo.wes.air.vol.T", newName="flo.wes.air.vol.T"),
      MatchVariable(name="flo.sou.air.vol.T", newName="flo.sou.air.vol.T"),
      MatchVariable(name="flo.nor.air.vol.T", newName="flo.nor.air.vol.T"),
      MatchVariable(name="flo.eas.air.vol.T", newName="flo.eas.air.vol.T"),
      MatchVariable(name="flo.cor.air.vol.T", newName="flo.cor.air.vol.T"),
      MatchVariable(name="flo.wes.air.senRelHum1.phi", newName="flo.wes.air.senRelHum.phi"),
      MatchVariable(name="flo.sou.air.senRelHum1.phi", newName="flo.sou.air.senRelHum.phi"),
      MatchVariable(name="flo.nor.air.senRelHum1.phi", newName="flo.nor.air.senRelHum.phi"),
      MatchVariable(name="flo.eas.air.senRelHum1.phi", newName="flo.eas.air.senRelHum.phi"),
      MatchVariable(name="flo.cor.air.senRelHum1.phi", newName="flo.cor.air.senRelHum.phi"),
      MatchVariable(name="flo.wes.air.vol.C[1]", newName="flo.wes.air.vol.C[1]"),
      MatchVariable(name="flo.sou.air.vol.C[1]", newName="flo.sou.air.vol.C[1]"),
      MatchVariable(name="flo.nor.air.vol.C[1]", newName="flo.nor.air.vol.C[1]"),
      MatchVariable(name="flo.eas.air.vol.C[1]", newName="flo.eas.air.vol.C[1]"),
      MatchVariable(name="flo.cor.air.vol.C[1]", newName="flo.cor.air.vol.C[1]"),
      MatchVariable(name="flo.wes.air.vol.C[2]", newName="flo.wes.air.vol.C[2]"),
      MatchVariable(name="flo.sou.air.vol.C[2]", newName="flo.sou.air.vol.C[2]"),
      MatchVariable(name="flo.nor.air.vol.C[2]", newName="flo.nor.air.vol.C[2]"),
      MatchVariable(name="flo.eas.air.vol.C[2]", newName="flo.eas.air.vol.C[2]"),
      MatchVariable(name="flo.cor.air.traceSubstance.C", newName="flo.cor.air.traceSubstance.C"),
      MatchVariable(name="flo.wes.air.traceSubstance.C", newName="flo.wes.air.traceSubstance.C"),
      MatchVariable(name="flo.sou.air.traceSubstance.C", newName="flo.sou.air.traceSubstance.C"),
      MatchVariable(name="flo.nor.air.traceSubstance.C", newName="flo.nor.air.traceSubstance.C"),
      MatchVariable(name="flo.eas.air.traceSubstance.C", newName="flo.eas.air.traceSubstance.C"),
      MatchVariable(name="flo.cor.air.vol.C[2]", newName="flo.cor.air.vol.C[2]"),
      MatchVariable(name="hvac.cor.terHea.Q1_flow", newName="hvac.cor.terHea.Q1_flow"),
      MatchVariable(name="hvac.nor.terHea.Q1_flow", newName="hvac.nor.terHea.Q1_flow"),
      MatchVariable(name="hvac.eas.terHea.Q1_flow", newName="hvac.eas.terHea.Q1_flow"),
      MatchVariable(name="hvac.wes.terHea.Q1_flow", newName="hvac.wes.terHea.Q1_flow"),
      MatchVariable(name="hvac.sou.terHea.Q1_flow", newName="hvac.sou.terHea.Q1_flow"),
      MatchVariable(name="flo.intGaiFra.y[1]", newName="flo.intGaiFra.y[1]"),
      MatchVariable(name="hvac.dpRetDuc.port_b.C_outflow[2]", newName="hvac.dpRetDuc.port_b.C_outflow[2]"),
      MatchVariable(name="hvac.dpRetDuc.port_b.C_outflow[1]", newName="hvac.dpRetDuc.port_b.C_outflow[1]"),
      MatchVariable(name="hvac.senSupFlo.port_b.C_outflow[2]", newName="hvac.senSupFlo.port_b.C_outflow[2]"),
      MatchVariable(name="hvac.senSupFlo.port_b.C_outflow[1]", newName="hvac.senSupFlo.port_b.C_outflow[1]"),
      MatchVariable(name="hvac.eco.port_Out.m_flow", newName="hvac.eco.port_Out.m_flow"),
      MatchVariable(name="hvac.eco.port_Sup.m_flow", newName="hvac.eco.port_Sup.m_flow"),
      MatchVariable(name="hvac.dpRetDuc.m_flow", newName="hvac.dpRetDuc.m_flow"),
      MatchVariable(name="flo.eas.ports[1].m_flow", newName="flo.eas.ports[1].m_flow"),
      MatchVariable(name="flo.wes.ports[1].m_flow", newName="flo.wes.ports[1].m_flow"),
      MatchVariable(name="flo.cor.ports[1].m_flow", newName="flo.cor.ports[1].m_flow"),
      MatchVariable(name="flo.nor.ports[1].m_flow", newName="flo.nor.ports[1].m_flow"),
      MatchVariable(name="flo.sou.ports[1].m_flow", newName="flo.sou.ports[1].m_flow"),
      MatchVariable(name="hvac.conFanSup.u_m", newName="hvac.conFanSup.u_m"),
      MatchVariable(name="hvac.conFanSup.y", newName="hvac.conFanSup.y"),
      MatchVariable(name="hvac.fanSup.port_a.p", newName="hvac.fanSup.port_a.p"),
      MatchVariable(name="hvac.TRet.T", newName="hvac.TRet.T"),
      MatchVariable(name="hvac.TSup.T", newName="hvac.TSup.T"),
      MatchVariable(name="weaBus.TDryBul", newName="weaBus.TDryBul"),
      MatchVariable(name="hvac.conEco.VOut_flow_min", newName="hvac.conEco.VOut_flow_min"),
      MatchVariable(name="hvac.pSetDuc.limPID.u_m", newName="hvac.pSetDuc.limPID.u_m"),
      MatchVariable(name="hvac.dpDisSupFan.p_rel", newName="hvac.dpDisSupFan.p_rel"),
      MatchVariable(name="hvac.controlBus.occupied", newName="hvac.controlBus.occupied"),
      MatchVariable(name="flo.sickPerson.y[1]", newName="flo.sickPerson.y[1]"),
      MatchVariable(name="hvac.fanSup.dp", newName="hvac.fanSup.dp"),
      MatchVariable(name="hvac.fanSup.VMachine_flow", newName="hvac.fanSup.VMachine_flow"),
      MatchVariable(name="hvac.fanSup.P", newName="hvac.fanSup.P"),
      MatchVariable(name="hvac.filt.dp_nominal", newName="hvac.filt.dp_nominal"),
      MatchVariable(name="flo.sou.air.gUV.yP_GUV", newName="flo.sou.air.gUV.yP_GUV"),
      MatchVariable(name="flo.sou.air.pAC.yP_PAC", newName="flo.sou.air.pAC.yP_PAC"),
      MatchVariable(name="flo.cor.air.gUV.yP_GUV", newName="flo.cor.air.gUV.yP_GUV"),
      MatchVariable(name="flo.cor.air.pAC.yP_PAC", newName="flo.cor.air.pAC.yP_PAC"),
      MatchVariable(name="flo.nor.air.gUV.yP_GUV", newName="flo.nor.air.gUV.yP_GUV"),
      MatchVariable(name="flo.nor.air.pAC.yP_PAC", newName="flo.nor.air.pAC.yP_PAC"),
      MatchVariable(name="flo.eas.air.gUV.yP_GUV", newName="flo.eas.air.gUV.yP_GUV"),
      MatchVariable(name="flo.eas.air.pAC.yP_PAC", newName="flo.eas.air.pAC.yP_PAC"),
      MatchVariable(name="flo.wes.air.gUV.yP_GUV", newName="flo.wes.air.gUV.yP_GUV"),
      MatchVariable(name="flo.wes.air.pAC.yP_PAC", newName="flo.wes.air.pAC.yP_PAC"),
      MatchVariable(name="flo.intGaiFra.y[1]", newName="flo.intGaiFra.y[1]"),
      MatchVariable(name="hvac.VAVBox[1].VSup_flow", newName="hvac.VAVBox[1].VSup_flow"),
      MatchVariable(name="hvac.VAVBox[2].VSup_flow", newName="hvac.VAVBox[2].VSup_flow"),
      MatchVariable(name="hvac.VAVBox[3].VSup_flow", newName="hvac.VAVBox[3].VSup_flow"),
      MatchVariable(name="hvac.VAVBox[4].VSup_flow", newName="hvac.VAVBox[4].VSup_flow"),
      MatchVariable(name="hvac.VAVBox[5].VSup_flow", newName="hvac.VAVBox[5].VSup_flow"),
      MatchVariable(name="hvac.VOut1.V_flow", newName="hvac.VOut1.V_flow"),
      MatchVariable(name="hvac.senSupFlo.V_flow", newName="hvac.senSupFlo.V_flow"),
      MatchVariable(name="hvac.VAVBox[1].terHea.Q2_flow", newName="hvac.VAVBox[1].terHea.Q2_flow"),
      MatchVariable(name="hvac.VAVBox[2].terHea.Q2_flow", newName="hvac.VAVBox[2].terHea.Q2_flow"),
      MatchVariable(name="hvac.VAVBox[3].terHea.Q2_flow", newName="hvac.VAVBox[3].terHea.Q2_flow"),
      MatchVariable(name="hvac.VAVBox[4].terHea.Q2_flow", newName="hvac.VAVBox[4].terHea.Q2_flow"),
      MatchVariable(name="hvac.VAVBox[5].terHea.Q2_flow", newName="hvac.VAVBox[5].terHea.Q2_flow"),
      MatchVariable(name="hvac.inDucGUV.yP_GUV", newName="hvac.inDucGUV.yP_GUV"),
      MatchVariable(name="flo.sou.flowPAC", newName="flo.sou.flowPAC"),
      MatchVariable(name="flo.cor.flowPAC", newName="flo.cor.flowPAC"),
      MatchVariable(name="flo.nor.flowPAC", newName="flo.nor.flowPAC"),
      MatchVariable(name="flo.eas.flowPAC", newName="flo.eas.flowPAC"),
      MatchVariable(name="flo.wes.flowPAC", newName="flo.wes.flowPAC"),
      MatchVariable(name="hvac.inDucGUV.filt.eff1", newName="hvac.inDucGUV.filt.eff1"),
      MatchVariable(name="hvac.heaCoi.Q1_flow", newName="hvac.heaCoi.Q1_flow"),
      MatchVariable(name="hvac.cooCoi.Q1_flow", newName="hvac.cooCoi.Q1_flow")})},
                  Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end OutputVariableSelection;

    model ViralDecay

      parameter Real kdec(min=0) = 0.48
        "Decay rate of virus";

      parameter Real V(min=0) = 100
        "Room volume";
      Modelica.Blocks.Math.Gain gain(k=-kdec*1.2*V/3600)
        annotation (Placement(transformation(extent={{-6,-10},{14,10}})));
      Modelica.Blocks.Interfaces.RealInput u
        annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
      Modelica.Blocks.Interfaces.RealOutput y
        annotation (Placement(transformation(extent={{100,-10},{120,10}})));
    equation
      connect(u, gain.u)
        annotation (Line(points={{-120,0},{-8,0}}, color={0,0,127}));
      connect(gain.y, y)
        annotation (Line(points={{15,0},{110,0}}, color={0,0,127}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p>This block calculates the concentration of a virus after viral decay using the decay rate of the virus(<span style=\"font-family: Courier New;\">kdec</span>) and the volume of the room(<span style=\"font-family: Courier New;\">V</span>) using the following equation:</p>
</html>"));
    end ViralDecay;

    partial model PartialInDuctGUV "Partial model for an in duct GUV"
        extends Buildings.Fluid.Interfaces.PartialTwoPortInterface(
         show_T=false,
         dp(nominal=if dp_nominal_pos > Modelica.Constants.eps
              then dp_nominal_pos else 1),
         m_flow(
            nominal=if m_flow_nominal_pos > Modelica.Constants.eps
              then m_flow_nominal_pos else 1),
         final m_flow_small = 1E-4*abs(m_flow_nominal));

      constant Boolean homotopyInitialization = true "= true, use homotopy method"
        annotation(HideResult=true);

      parameter Boolean from_dp = false
        "= true, use m_flow = f(dp) else dp = f(m_flow)"
        annotation (Evaluate=true, Dialog(tab="Advanced"));

      parameter Modelica.Units.SI.PressureDifference dp_nominal(displayUnit="Pa")
        "Pressure drop at nominal mass flow rate"
        annotation (Dialog(group="Nominal condition"));

      parameter Boolean linearized = false
        "= true, use linear relation between m_flow and dp for any flow rate"
        annotation(Evaluate=true, Dialog(tab="Advanced"));

      parameter Modelica.Units.SI.MassFlowRate m_flow_turbulent(min=0)
        "Turbulent flow if |m_flow| >= m_flow_turbulent";

      parameter Real eff(min=0, max=1) = 0.85
        "Efficiency of HVAC filter";

    protected
      parameter Medium.ThermodynamicState sta_default=
         Medium.setState_pTX(T=Medium.T_default, p=Medium.p_default, X=Medium.X_default);
      parameter Modelica.Units.SI.DynamicViscosity eta_default=
          Medium.dynamicViscosity(sta_default)
        "Dynamic viscosity, used to compute transition to turbulent flow regime";

      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal_pos=abs(
          m_flow_nominal) "Absolute value of nominal flow rate";
      final parameter Modelica.Units.SI.PressureDifference dp_nominal_pos(
          displayUnit="Pa") = abs(dp_nominal)
        "Absolute value of nominal pressure difference";
      Delays.DelayFirstOrder                 vol(
        redeclare final package Medium = Medium,
        use_C_flow=false,
        final tau=1,
        final energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
        final m_flow_nominal=m_flow_nominal,
        final m_flow_small=m_flow_small,
        final prescribedHeatFlowRate=true,
        final allowFlowReversal=allowFlowReversal,
        nPorts=2) "Fluid volume for dynamic model"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    initial equation
      assert(homotopyInitialization, "In " + getInstanceName() +
        ": The constant homotopyInitialization has been modified from its default value. This constant will be removed in future releases.",
        level = AssertionLevel.warning);

    equation
      // Isenthalpic state transformation (no storage and no loss of energy)
      //port_a.h_outflow = if allowFlowReversal then inStream(port_b.h_outflow) else Medium.h_default;
      //port_b.h_outflow = inStream(port_a.h_outflow);

      // Mass balance (no storage)
      //port_a.m_flow + port_b.m_flow = 0;

      // Transport of substances
      //port_a.Xi_outflow = if allowFlowReversal then inStream(port_b.Xi_outflow) else Medium.X_default[1:Medium.nXi];
      //port_b.Xi_outflow = inStream(port_a.Xi_outflow);

      //port_a.C_outflow = if allowFlowReversal then inStream(port_b.C_outflow) else zeros(Medium.nC);
      //port_b.C_outflow[2] = (1-eff)*inStream(vol.ports[2].C_outflow[2]);
      //port_b.C_outflow[1] = inStream(vol.ports[2].C_outflow[1]);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Rectangle(
              extent=DynamicSelect({{-100,10},{-100,10}}, {{100,10},{100+200*max(-1, min(0, m_flow/(abs(m_flow_nominal)))),-10}}),
              lineColor={28,108,200},
              fillColor={255,0,0},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Rectangle(
              extent=DynamicSelect({{-100,10},{-100,10}}, {{-100,10},{-100+200*min(1, max(0, m_flow/abs(m_flow_nominal))),-10}}),
              lineColor={28,108,200},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None)}),
              defaultComponentName="res",
    Documentation(info="<html>
<p>
Partial model for a flow resistance, possible with variable flow coefficient.
Models that extend this class need to implement an equation that relates
<code>m_flow</code> and <code>dp</code>, and they need to assign the parameter
<code>m_flow_turbulent</code>.
</p>
<p>
See for example
<a href=\"modelica://Buildings.Fluid.FixedResistances.PressureDrop\">
Buildings.Fluid.FixedResistances.PressureDrop</a> for a model that extends
this base class.
</p>
</html>",     revisions="<html>
<ul>
<li>
April 14, 2020, by Michael Wetter:<br/>
Changed <code>homotopyInitialization</code> to a constant.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1341\">IBPSA, #1341</a>.
</li>
<li>
February 26, 2020, by Michael Wetter:<br/>
Changed icon to display its operating state.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1294\">#1294</a>.
</li>
<li>
October 25, 2019, by Jianjun Hu:<br/>
Improved icon graphics annotation. This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1225\">#1225</a>.
</li>
<li>
November 3, 2016, by Michael Wetter:<br/>
Removed start value for pressure difference
to simplify the parameter window.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/552\">#552</a>.
</li>
<li>
January 26, 2016, by Michael Wetter:<br/>
Avoided assignment of <code>dp(nominal=0)</code> if <code>dp_nominal_pos = 0</code>
and of <code>m_flow(nominal=0)</code> if <code>m_flow_nominal_pos = 0</code>
as nominal values are not allowed to be zero.
</li>
<li>
January 22, 2016, by Michael Wetter:<br/>
Corrected type declaration of pressure difference.
This is
for <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/404\">#404</a>.
</li>
<li>
August 15, 2015, by Filip Jorissen:<br/>
Implemented more efficient computation of <code>port_a.Xi_outflow</code>,
<code>port_a.h_outflow</code>
and <code>port_a.C_outflow</code> when <code>allowFlowReversal=false</code>.
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/281\">#281</a>.
</li>
<li>
January 13, 2015, by Marcus Fuchs:<br/>
Revised revisions section (there were two revisions statements)
</li>
<li>
November 20, 2014 by Michael Wetter:<br/>
Removed <code>start</code> attribute for <code>m_flow</code>
as this is already set in its base class.
</li>
<li>
October 8, 2013 by Michael Wetter:<br/>
Removed propagation of <code>show_V_flow</code>
to base class as it has no longer this parameter.
</li>
<li>
December 14, 2012 by Michael Wetter:<br/>
Renamed protected parameters for consistency with the naming conventions.
</li>
<li>
February 12, 2012, by Michael Wetter:<br/>
Removed duplicate declaration of <code>m_flow_nominal</code>.
</li>
<li>
February 3, 2012, by Michael Wetter:<br/>
Made assignment of <code>m_flow_small</code> <code>final</code> as it is no
longer used in the base class.
</li>
<li>
January 16, 2012, by Michael Wetter:<br/>
To simplify object inheritance tree, revised base classes
<code>Buildings.Fluid.BaseClasses.PartialResistance</code>,
<code>Buildings.Fluid.Actuators.BaseClasses.PartialTwoWayValve</code>,
<code>Buildings.Fluid.Actuators.BaseClasses.PartialDamperExponential</code>,
<code>Buildings.Fluid.Actuators.BaseClasses.PartialActuator</code>
and model
<code>Buildings.Fluid.FixedResistances.PressureDrop</code>.
</li>
<li>
August 5, 2011, by Michael Wetter:<br/>
Moved linearized pressure drop equation from the function body to the equation
section. With the previous implementation,
the symbolic processor may not rearrange the equations, which can lead
to coupled equations instead of an explicit solution.
</li>
<li>
June 20, 2011, by Michael Wetter:<br/>
Set start values for <code>m_flow</code> and <code>dp</code> to zero, since
most HVAC systems start at zero flow. With this change, the start values
appear in the GUI and can be set by the user.
</li>
<li>
April 2, 2011 by Michael Wetter:<br/>
Added <code>m_flow_nominal_pos</code> and <code>dp_nominal_pos</code> to allow
providing negative nominal values which will be used, for example, to set start
values of flow splitters which may have negative flow rates and pressure drop
at the initial condition.
</li>
<li>
March 27, 2011, by Michael Wetter:<br/>
Added <code>homotopy</code> operator.
</li>
<li>
March 23, 2011 by Michael Wetter:<br/>
Added homotopy operator.
</li>
<li>
March 30, 2010 by Michael Wetter:<br/>
Changed base classes to allow easier initialization.
</li>
<li>
April 13, 2009, by Michael Wetter:<br/>
Extracted pressure drop computation and implemented it in the
new model
<a href=\"modelica://Buildings.Fluid.BaseClasses.FlowModels.BasicFlowModel\">
Buildings.Fluid.BaseClasses.FlowModels.BasicFlowModel</a>.
</li>
<li>
September 18, 2008, by Michael Wetter:<br/>
Added equations for the mass balance of extra species flow,
i.e., <code>C</code> and <code>mC_flow</code>.
</li>
<li>
July 20, 2007 by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"));
    end PartialInDuctGUV;
  end BaseClasses;

  package Examples
    model inDucGUVExample
      replaceable package MediumA = Buildings.Media.Air(extraPropertiesNames={"CO2", "COVID"}) "Medium model for air";
      parameter Modelica.Units.SI.MassFlowRate mAir_flow_nominal = 6.71
        "Design mass flow rate";
      Buildings.Fluid.FixedResistances.InDuctGUV
                                       inDucGUV(
        dp_nominal=0,
        eff=1,
        kpow=200,
        m_flow_nominal=mAir_flow_nominal,
        redeclare package Medium = MediumA)
        annotation (Placement(transformation(extent={{-10,-12},{10,8}})));
      Buildings.Fluid.Sources.MassFlowSource_T boundary(
        use_C_in=false,
        C={0,1},
        use_m_flow_in=true,
        redeclare package Medium = MediumA,
        nPorts=1) annotation (Placement(transformation(extent={{-58,-12},{-38,8}})));
      Buildings.Fluid.Sources.Boundary_pT bou(nPorts=1, redeclare package
          Medium =                                                                 MediumA)
        annotation (Placement(transformation(extent={{72,-10},{52,10}})));
      Modelica.Blocks.Sources.Ramp ramp(
        height=0.75*mAir_flow_nominal,
        duration=3600,
        offset=0.25*mAir_flow_nominal,
        startTime=120)
        annotation (Placement(transformation(extent={{-88,26},{-68,46}})));
      Modelica.Blocks.Sources.BooleanConstant booleanConstant
        annotation (Placement(transformation(extent={{-66,-46},{-46,-26}})));
    equation
      connect(boundary.ports[1], inDucGUV.port_a)
        annotation (Line(points={{-38,-2},{-10,-2}}, color={0,127,255}));
      connect(inDucGUV.port_b, bou.ports[1]) annotation (Line(points={{10,-2},{46,-2},
              {46,0},{52,0}}, color={0,127,255}));
      connect(booleanConstant.y, inDucGUV.u) annotation (Line(points={{-45,-36},{-20,
              -36},{-20,-10},{-12,-10}}, color={255,0,255}));
      connect(ramp.y, boundary.m_flow_in) annotation (Line(points={{-67,36},{-62,36},
              {-62,12},{-66,12},{-66,6},{-60,6}}, color={0,0,127}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)));
    end inDucGUVExample;

    model PACExample
      import Buildings.AirCleaning;
      replaceable package Medium = Buildings.Media.Air(extraPropertiesNames={"CO2", "COVID"}) "Medium model for air";
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6.71
        "Design mass flow rate";
      Buildings.Fluid.MixingVolumes.MixingVolumeMoistAir
                                               vol(
        redeclare package Medium = Medium,
        final energyDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial,
        final massDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial,
        final V=500,
        final C_start={0,1},
        final mSenFac=1,
        final m_flow_nominal=m_flow_nominal,
        final prescribedHeatFlowRate=true,
        final nPorts=2,
        m_flow_small=1E-4*abs(m_flow_nominal),
        allowFlowReversal=true,
        final use_C_flow=true)        "Room air volume"
        annotation (Placement(transformation(extent={{8,-8},{28,12}})));
      Modelica.Blocks.Sources.Constant conZero(k=0)
        annotation (Placement(transformation(extent={{-60,22},{-40,42}})));
      AirCleaning.PAC pAC(
        redeclare package Medium = Medium,
        eff={0,0.9997},
        nPACs=2,
        flowPAC=400*1.2/3600,
        kpow=100)
        annotation (Placement(transformation(extent={{-36,-62},{-16,-42}})));
      Modelica.Blocks.Sources.RealExpression patConc(y=vol.C[2])
        "virus concentration"
        annotation (Placement(transformation(extent={{-90,-40},{-70,-20}})));
      Buildings.Fluid.Sources.MassFlowSource_T sup(
        use_C_in=false,
        C={0,1},
        use_m_flow_in=false,
        redeclare package Medium = Medium,
        m_flow=0*500*1.2/3600,
        nPorts=1)
        annotation (Placement(transformation(extent={{-68,-10},{-48,10}})));
      Buildings.Fluid.Sources.Boundary_pT bou(nPorts=1, redeclare package
          Medium =
            Medium)
        annotation (Placement(transformation(extent={{58,-38},{38,-18}})));
      Modelica.Blocks.Sources.BooleanConstant ConPACEna "True when PAC is on"
        annotation (Placement(transformation(extent={{-90,-72},{-70,-52}})));
      AirCleaning.BaseClasses.ViralDecay viralDecay
        annotation (Placement(transformation(extent={{-36,-90},{-16,-70}})));
      Modelica.Blocks.Sources.RealExpression patConc1(y=vol.C[1])
        "CO2 concentration"
        annotation (Placement(transformation(extent={{-62,-30},{-42,-10}})));
    equation
      connect(conZero.y, vol.mWat_flow) annotation (Line(points={{-39,32},{-22,32},
              {-22,10},{6,10}}, color={0,0,127}));
      connect(pAC.diss, vol.heatPort) annotation (Line(points={{-15.6,-60},{-6,-60},
              {-6,2},{8,2}},
                    color={191,0,0}));
      connect(sup.ports[1], vol.ports[1]) annotation (Line(points={{-48,0},{-4,0},{
              -4,-12},{16,-12},{16,-8}}, color={0,127,255}));
      connect(vol.ports[2], bou.ports[1])
        annotation (Line(points={{20,-8},{20,-28},{38,-28}},  color={0,127,255}));
      connect(ConPACEna.y, pAC.uPACEna) annotation (Line(points={{-69,-62},{-46,-62},
              {-46,-54},{-38,-54}}, color={255,0,255}));
      connect(patConc.y, viralDecay.u) annotation (Line(points={{-69,-30},{-44,-30},
              {-44,-80},{-38,-80}}, color={0,0,127}));
      connect(pAC.yC_flow[1], vol.C_flow[1]) annotation (Line(points={{-15,-48},{0,
              -48},{0,-10},{2,-10},{2,-4},{6,-4}}, color={0,0,127}));
      connect(pAC.yC_flow[2], vol.C_flow[2]) annotation (Line(points={{-15,-48},{0,
              -48},{0,-10},{2,-10},{2,-4},{6,-4}}, color={0,0,127}));
      connect(patConc1.y, pAC.C[1]) annotation (Line(points={{-41,-20},{-38,-20},{
              -38,-38},{-46,-38},{-46,-48},{-38,-48}}, color={0,0,127}));
      connect(patConc.y, pAC.C[2]) annotation (Line(points={{-69,-30},{-44,-30},{
              -44,-48},{-38,-48}}, color={0,0,127}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p>This example simulates <a href=\"modelica://Buildings.AirCleaning.PAC\">AirCleaning.PAC</a> for the PAC outflow concentration of a virus and CO<sub>2. </sub></p>
</html>"));
    end PACExample;
  end Examples;
  annotation ();
end AirCleaning;
