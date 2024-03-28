within Buildings;
package AirCleaning

  model PAC "In-room portable air cleaner"

    replaceable package Medium =
      Modelica.Media.Interfaces.PartialMedium "Medium in the component";

    parameter Real eff[Medium.nC](min=0, max=1) = 0.9997
      "trace species removal efficiency";

    parameter Modelica.Units.SI.MassFlowRate flow_PAC(min=0) = 0.094
      "PAC flow rate";

    parameter Modelica.Units.SI.Power kpow(min=0) = 50
      "Rated power";

    Modelica.Blocks.Interfaces.RealInput C[Medium.nC] "Zone concentration(s)"
      annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
    Modelica.Blocks.Interfaces.RealOutput yC_flow[Medium.nC] "Concentration outflow"
      annotation (Placement(transformation(extent={{100,30},{120,50}})));
    Modelica.Blocks.Interfaces.BooleanInput uPACEna "True when PAC is on"
      annotation (Placement(transformation(extent={{-140,-40},{-100,0}})));
    Modelica.Blocks.Math.BooleanToReal booleanToReal
      annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a diss "Heat dissipation"
      annotation (Placement(transformation(extent={{94,-90},{114,-70}})));

  protected
    Modelica.Blocks.Math.Gain pow(final k=kpow)       "power of PAC"
      annotation (Placement(transformation(extent={{-20,-30},{0,-10}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prePow(final alpha=0)
      "Prescribed power of the PAC dissipated as heat"
      annotation (Placement(transformation(extent={{38,-30},{58,-10}})));

  equation
    connect(uPACEna, booleanToReal.u)
      annotation (Line(points={{-120,-20},{-82,-20}}, color={255,0,255}));
    connect(booleanToReal.y, pow.u)
      annotation (Line(points={{-59,-20},{-22,-20}}, color={0,0,127}));
    connect(pow.y, prePow.Q_flow) annotation (Line(points={{1,-20},{38,-20}},
                                       color={0,0,127}));
    connect(prePow.port, diss) annotation (Line(points={{58,-20},{90,-20},{90,-80},
            {104,-80}}, color={191,0,0}));

    // calculate trace species mass flow exiting PAC
    for i in 1:Medium.nC loop
      yC_flow[i]  = booleanToReal.y*(-flow_PAC)*eff[i]*C[i];
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
<p>This model is designed to simulate the removal of trace species using a portable air cleaner (PAC). The PAC model calculates the mass removal rate of trace species by the PAC (<span style=\"font-family: Courier New;\">yC_flow[]</span>) based on the equation described below when the PAC is on (<span style=\"font-family: Courier New;\">uPACEna</span>). The PACs are modeled to consume energy using a constant power rating (<span style=\"font-family: Courier New;\">kpow</span>) and heat is dissipated into the zone based on the consumed power. </p>
<h4>Main Equations </h4>
<p align=\"center\"><span style=\"font-family: Courier New;\">Ċ</span><sub>PAC</sub> = eff<sub>PAC</sub>* <span style=\"font-family: Courier New;\">flow</span><sub>PAC</sub>*c<sub>zone</sub> </p>
<p>where <span style=\"font-family: Courier New;\">Ċ</span><sub>PAC</sub> is the rate of trace species removal by the PAC, eff<sub>PAC</sub> is the trace species removal efficency, <span style=\"font-family: Courier New;\">flow</span><sub>PAC</sub> is the mass airflow rate of the PAC, and c<sub>zone</sub> is the trace species concentration in the zone where the PAC is located.</p>
<h4>Assumptions</h4>
<p>The model assumes well-mixed zones with uniform concentrations. </p>
</html>"));
  end PAC;

  model RoomGUV "In-room GUV"

    replaceable package Medium =
      Modelica.Media.Interfaces.PartialMedium "Medium in the component";

    parameter Real frad(min=0, max=1) = 0.2
      "Fraction of irradiated space";

    parameter Modelica.Units.SI.Irradiance Eavg=0.5
      "Effluence rate";

    parameter Real krad[Medium.nC](unit="m2/J")=0.28
      "Inactivation constant";

    parameter Modelica.Units.SI.Power kpow(min=0) = 120
      "Rated power";

    parameter Modelica.Units.SI.Volume V=120
      "Zone volume";

    Modelica.Blocks.Interfaces.RealInput C[Medium.nC] "Zone concentration"
      annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
    Modelica.Blocks.Interfaces.RealOutput yC_flow[Medium.nC]
      "Concentration outflow"
      annotation (Placement(transformation(extent={{100,30},{120,50}})));
    Modelica.Blocks.Interfaces.BooleanInput u "on/off"
      annotation (Placement(transformation(extent={{-140,-40},{-100,0}})));
    Modelica.Blocks.Math.BooleanToReal booleanToReal
      annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a diss "Heat dissipation"
      annotation (Placement(transformation(extent={{94,-90},{114,-70}})));

  protected
    Modelica.Blocks.Math.Gain pow(final k=kpow)
                                               "power of GUV"
      annotation (Placement(transformation(extent={{-20,-30},{0,-10}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prePow(final alpha=0)
      "Prescribed power (=heat and flow work) flow for dynamic model"
      annotation (Placement(transformation(extent={{42,-30},{62,-10}})));

  equation

    for i in 1:Medium.nC loop
      yC_flow[i]  = booleanToReal.y*1.2*(-frad)*Eavg*krad[i]*V*C[i];
    end for;

    connect(u, booleanToReal.u)
      annotation (Line(points={{-120,-20},{-82,-20}}, color={255,0,255}));
    connect(booleanToReal.y, pow.u)
      annotation (Line(points={{-59,-20},{-22,-20}}, color={0,0,127}));
    connect(pow.y, prePow.Q_flow) annotation (Line(points={{1,-20},{42,-20}},
                                       color={0,0,127}));
    connect(prePow.port, diss) annotation (Line(points={{62,-20},{90,-20},{90,-80},
            {104,-80}}, color={191,0,0}));
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
<p>This model is designed to simulate the inactivation of trace species by either upper- or whole-room germicidal ultraviolet (GUV) devices when the GUV is enabled <span style=\"font-family: Courier New;\">uGUVEna</span>. The RoomGUV model calculates the mass inactivation rate of trace species in the zone from the GUV(yC_flow) based on the equation described below. The GUV device consumes energy based on a constant power rating, <span style=\"font-family: Courier New;\">kpow</span>, which is ultimately dissipated into the zone as heat. </p>
<h4>Main Equations </h4>
<p>Ċ<sub>GUV</sub> = E<sub>avg</sub> k<sub>rad</sub>f<sub>rad</sub>V<sub>zone</sub>c<sub>zone</sub> </p>
<p>where Ċ<sub>PAC</sub> is the inactivation rate of trace species by the GUV, E<sub>avg</sub> is the average fluence rate of the GUV device, k<sub>rad</sub> is susceptibility of the trace species to the GUV irrdiation, f<sub>rad</sub> is the fraction of irradiated volume in the zone, V<sub>zone</sub> is the total volume of the zone, and c<sub>zone</sub> is the trace species concentration in the zone where the GUV is located. </p>
<h4>Assumptions </h4>
<p>This assumes that the zone is well mixed and the parameters E<sub>avg</sub>,k<sub>rad</sub>, and f<sub>rad</sub> are constant for the given simulation when the GUV device is on. </p>
</html>"));
  end RoomGUV;

  model DuctGUV "In Duct GUV"
    extends Buildings.AirCleaning.BaseClasses.PartialDuctGUV(
      final m_flow_turbulent = if computeFlowResistance then deltaM * m_flow_nominal_pos else 0, vol(
          nPorts=2));

    parameter Real deltaM(min=1E-6) = 0.3
      "Fraction of nominal mass flow rate where transition to turbulent occurs"
         annotation(Evaluate=true,
                    Dialog(group = "Transition to laminar",
                           enable = not linearized));

    parameter Real kpow(min=0) = 120
      "Rated power";

    parameter Real kGUV[Medium.nC](min=0)
      "Inactivation constant";

    parameter Boolean addPowerToMedium=true
      "Set to false to avoid any power (=heat and flow work) being added to medium (may give simpler equations)";

    final parameter Real k = if computeFlowResistance then
          m_flow_nominal_pos / sqrt(dp_nominal_pos) else 0
      "Flow coefficient, k=m_flow/sqrt(dp), with unit=(kg.m)^(1/2)";
    Buildings.AirCleaning.BaseClasses.InDuctGUVCalc guvCal(
      redeclare package Medium = Medium,
      m_flow_nominal=m_flow_nominal,
      dp_nominal=dp_nominal,
      kGUV=kGUV) annotation (Placement(transformation(extent={{44,-10},{64,10}})));
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
    Modelica.Blocks.Math.Gain pGUV(final k=kpow) "power of GUV"
      annotation (Placement(transformation(extent={{-48,-60},{-28,-40}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prePow(final alpha=0)
   if addPowerToMedium
      "Prescribed power (=heat and flow work) flow for dynamic model"
      annotation (Placement(transformation(extent={{-20,-60},{0,-40}})));
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

    connect(pGUV.y, prePow.Q_flow)
      annotation (Line(points={{-27,-50},{-20,-50}}, color={0,0,127}));
    connect(guvCal.port_b, port_b)
      annotation (Line(points={{64,0},{100,0}}, color={0,127,255}));
    connect(booleanToReal.y, pGUV.u) annotation (Line(points={{-59,-80},{-54,-80},
            {-54,-50},{-50,-50}}, color={0,0,127}));
    connect(u, guvCal.u) annotation (Line(points={{-120,-80},{-90,-80},{-90,-18},
            {14,-18},{14,-8},{42,-8}}, color={255,0,255}));
    connect(prePow.port, vol.heatPort) annotation (Line(points={{0,-50},{0,0}},
                                   color={191,0,0}));
    connect(port_a, vol.ports[1]) annotation (Line(points={{-100,0},{-6,0},{-6,
            -12},{2,-12},{2,-14},{8,-14},{8,-10}}, color={0,127,255}));
    connect(vol.ports[2], guvCal.port_a) annotation (Line(points={{12,-10},{56,
            -10},{56,0},{44,0}}, color={0,127,255}));
    annotation (defaultComponentName="res",
  Documentation(info="<html>
<p>Model of an in-duct GUV. </p>
<h4>Assumptions</h4>
<h4>Important parameters</h4>
<h4>Notes</h4>
<h4>Implementation</h4>
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
  end DuctGUV;

  package BaseClasses

    partial model PartialDuctGUV "Partial model for an in duct GUV"
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

      parameter Real kGUV[Medium.nC](min=0) = 1
        "Inactivation constant";

      Modelica.Blocks.Math.BooleanToReal booleanToReal
        annotation (Placement(transformation(extent={{-80,-90},{-60,-70}})));
      Modelica.Blocks.Interfaces.BooleanInput u "on/off"
        annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));
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
      Buildings.Fluid.Delays.DelayFirstOrder                 vol(
        redeclare final package Medium = Medium,
        use_C_flow=false,
        final tau=1,
        final energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
        final m_flow_nominal=m_flow_nominal,
        final m_flow_small=m_flow_small,
        final prescribedHeatFlowRate=true,
        final allowFlowReversal=allowFlowReversal,
        nPorts=2) "Fluid volume for dynamic model"
        annotation (Placement(transformation(extent={{0,-10},{20,10}})));
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
      //for i in 1:Medium.nC loop
       // port_b.C_outflow[i] = booleanToReal.y*(1 - exp(-kGUV[i]*m_flow_nominal/m_flow))*inStream(vol.ports[2].C_outflow[i]);
      //end for;
      //port_b.C_outflow[2] = inStream(port_a.C_outflow[2]);
      //port_b.C_outflow[1] = inStream(port_a.C_outflow[1]);

      connect(u,booleanToReal. u)
        annotation (Line(points={{-120,-80},{-82,-80}}, color={255,0,255}));
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
<p>Partial model for duct GUV.</p>
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
    end PartialDuctGUV;

    partial model PartialInDuctGUVCalc "Partial model for an in duct GUV"
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

      parameter Modelica.Units.SI.MassFlowRate m_flow_turbulent
        "Turbulent flow if |m_flow| >= m_flow_turbulent";

      parameter Real kGUV[Medium.nC](min=0)
        "Inactivation constant";

      Modelica.Units.SI.MassFlowRate m_flow_GUV(min = 0)
        "comment";

      Modelica.Blocks.Math.BooleanToReal booleanToReal
        annotation (Placement(transformation(extent={{-80,-90},{-60,-70}})));
      Modelica.Blocks.Interfaces.BooleanInput u "on/off"
        annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));
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
    initial equation
      assert(homotopyInitialization, "In " + getInstanceName() +
        ": The constant homotopyInitialization has been modified from its default value. This constant will be removed in future releases.",
        level = AssertionLevel.warning);

    equation
      // Isenthalpic state transformation (no storage and no loss of energy)
      port_a.h_outflow = if allowFlowReversal then inStream(port_b.h_outflow) else Medium.h_default;
      port_b.h_outflow = inStream(port_a.h_outflow);

      // Mass balance (no storage)
      port_a.m_flow + port_b.m_flow = 0;

      // Transport of substances
      port_a.Xi_outflow = if allowFlowReversal then inStream(port_b.Xi_outflow) else Medium.X_default[1:Medium.nXi];
      port_b.Xi_outflow = inStream(port_a.Xi_outflow);

      port_a.C_outflow = if allowFlowReversal then inStream(port_b.C_outflow) else zeros(Medium.nC);
      if m_flow<m_flow_small then
        m_flow_GUV = m_flow_small;
      else
        m_flow_GUV = m_flow;
      end if;
      for i in 1:Medium.nC loop
        port_b.C_outflow[i] = booleanToReal.y*(1 - exp(-kGUV[i]*m_flow_nominal/m_flow_GUV))*inStream(port_a.C_outflow[i]);
      end for;
      //port_b.C_outflow[2] = (1-eff)*inStream(vol.ports[2].C_outflow[2]);
      //port_b.C_outflow[1] = inStream(vol.ports[2].C_outflow[1]);

      connect(u,booleanToReal. u)
        annotation (Line(points={{-120,-80},{-82,-80}}, color={255,0,255}));
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
Partial model for a duct GUV trace species inactivation calculation.
</p>
<p>
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
    end PartialInDuctGUVCalc;

    model InDuctGUVCalc "HVAC filter"
      extends Buildings.AirCleaning.BaseClasses.PartialInDuctGUVCalc(
        final m_flow_turbulent = if computeFlowResistance then deltaM * m_flow_nominal_pos else 0);

      parameter Real deltaM(min=1E-6) = 0.3
        "Fraction of nominal mass flow rate where transition to turbulent occurs"
           annotation(Evaluate=true,
                      Dialog(group = "Transition to laminar",
                             enable = not linearized));
      final parameter Real k = if computeFlowResistance then
            m_flow_nominal_pos / sqrt(dp_nominal_pos) else 0
        "Flow coefficient, k=m_flow/sqrt(dp), with unit=(kg.m)^(1/2)";
    protected
      final parameter Boolean computeFlowResistance=(dp_nominal_pos > Modelica.Constants.eps)
        "Flag to enable/disable computation of flow resistance"
       annotation(Evaluate=true);
      final parameter Real coeff=
        if linearized and computeFlowResistance
        then if from_dp then k^2/m_flow_nominal_pos else m_flow_nominal_pos/k^2
        else 0
        "Precomputed coefficient to avoid division by parameter";
    initial equation
     if computeFlowResistance then
       assert(m_flow_turbulent > 0, "m_flow_turbulent must be bigger than zero.");
     end if;

     assert(m_flow_nominal_pos > 0, "m_flow_nominal_pos must be non-zero. Check parameters.");
    equation
      // Pressure drop calculation
      if computeFlowResistance then
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
      end if;  // computeFlowResistance

      annotation (defaultComponentName="res",
    Documentation(info="<html>
<p>Duct GUV calculation.</p>
<p><b>Assumptions</b> </p>
<h4>Important parameters</h4>
<h4>Notes</h4>
<h4>Implementation</h4>
</html>",     revisions="<html>
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
              extent={{-74,80},{72,62}},
              lineColor={28,108,200},
              fillColor={85,170,255},
              fillPattern=FillPattern.Solid),
            Line(points={{-60,56},{-60,-28}}, color={28,108,200}),
            Line(points={{-40,56},{-40,-28}}, color={28,108,200}),
            Line(points={{0,56},{0,-28}}, color={28,108,200}),
            Line(points={{60,56},{60,-28}}, color={28,108,200}),
            Line(points={{40,56},{40,-28}}, color={28,108,200}),
            Line(points={{20,56},{20,-28}}, color={28,108,200}),
            Line(points={{-20,56},{-20,-28}}, color={28,108,200})}));
    end InDuctGUVCalc;
  end BaseClasses;

  package Examples

    model PACExample
      import Buildings.AirCleaning;
      replaceable package Medium = Buildings.Media.Air(extraPropertiesNames={"CO2", "SARS-CoV-2"}) "Medium model for air";
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal=1.2*vol.V/3600
        "Design mass flow rate";
      Buildings.Fluid.MixingVolumes.MixingVolumeMoistAir
                                               vol(
        redeclare package Medium = Medium,
        final energyDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial,
        final massDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial,
        final V=500,
        final C_start={0,0},
        final mSenFac=1,
        final m_flow_nominal=m_flow_nominal,
        final prescribedHeatFlowRate=true,
        final nPorts=1,
        m_flow_small=1E-4*abs(m_flow_nominal),
        allowFlowReversal=true,
        final use_C_flow=true)        "Room air volume"
        annotation (Placement(transformation(extent={{40,-10},{60,10}})));
      Modelica.Blocks.Sources.Constant conZero(k=0)
        annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
      AirCleaning.PAC pAC(
        redeclare package Medium = Medium,
        eff={0,0.99},
        flow_PAC=100*1.2/3600,
        kpow=50)
        annotation (Placement(transformation(extent={{-30,-60},{-10,-40}})));
      Modelica.Blocks.Sources.RealExpression speConc2(y=vol.C[2])
        "SARS-CoV-2 concentration"
        annotation (Placement(transformation(extent={{-90,-50},{-70,-30}})));
      Buildings.Fluid.Sources.MassFlowSource_T sup(
        use_C_in=false,
        C={0.1,0.1},
        use_m_flow_in=false,
        redeclare package Medium = Medium,
        m_flow=m_flow_nominal,
        nPorts=1) "Supply airflow"
        annotation (Placement(transformation(extent={{-68,-10},{-48,10}})));
      Modelica.Blocks.Sources.BooleanConstant conPACEna "True when PAC is on"
        annotation (Placement(transformation(extent={{-90,-72},{-70,-52}})));
      Modelica.Blocks.Sources.RealExpression speConc1(y=vol.C[1])
        "CO2 concentration"
        annotation (Placement(transformation(extent={{-90,-30},{-70,-10}})));
    equation
      connect(pAC.diss, vol.heatPort) annotation (Line(points={{-9.6,-58},{-2,-58},{
              -2,2},{34,2},{34,0},{40,0}},
                    color={191,0,0}));
      connect(sup.ports[1], vol.ports[1]) annotation (Line(points={{-48,0},{6,0},{6,
              -16},{50,-16},{50,-10}},   color={0,127,255}));
      connect(conPACEna.y, pAC.uPACEna) annotation (Line(points={{-69,-62},{-46,-62},
              {-46,-52},{-32,-52}}, color={255,0,255}));
      connect(pAC.yC_flow[1], vol.C_flow[1]) annotation (Line(points={{-9,-46},{0,-46},
              {0,-6},{38,-6}},                     color={0,0,127}));
      connect(pAC.yC_flow[2], vol.C_flow[2]) annotation (Line(points={{-9,-46},{8,-46},
              {8,-6},{38,-6}},                     color={0,0,127}));
      connect(speConc1.y, pAC.C[1]) annotation (Line(points={{-69,-20},{-46,-20},{-46,
              -46},{-32,-46}},                         color={0,0,127}));
      connect(speConc2.y, pAC.C[2]) annotation (Line(points={{-69,-40},{-46,-40},{-46,
              -46},{-32,-46}}, color={0,0,127}));
      connect(conZero.y, vol.mWat_flow)
        annotation (Line(points={{-19,50},{38,50},{38,8}}, color={0,0,127}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p>This example simulates <a href=\"modelica://Buildings.AirCleaning.PAC\">AirCleaning.PAC</a> for a scenario with CO<sub>2. </sub>and SARS-CoV-2 as trace species. Airflow with a constant rate and concentration of trace species is supplied to a mixing volume. A PAC is connected to the mixing volume to simulate removal of trace species from the volume, where the PAC cannot remove the gaseous CO<sub>2</sub> but filters SARS-CoV-2 with an efficiency of 99&percnt;.</p>
</html>"),
        experiment(StopTime=3600, __Dymola_Algorithm="Dassl"));
    end PACExample;

    model DuctGUVExample
      replaceable package Medium = Buildings.Media.Air(extraPropertiesNames={"CO2", "SARS-CoV-2"}) "Medium model for air";
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6.71
        "Design mass flow rate";
      Buildings.AirCleaning.DuctGUV  ducGUV(
        dp_nominal=0,
        kGUV={10e6,0.69},
        kpow=200,
        m_flow_nominal=m_flow_nominal,
        redeclare package Medium = Medium)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Buildings.Fluid.Sources.MassFlowSource_T flo(
        use_C_in=false,
        C={1,1},
        use_m_flow_in=true,
        redeclare package Medium = Medium,
        nPorts=1)
        annotation (Placement(transformation(extent={{-48,-10},{-28,10}})));
      Modelica.Blocks.Sources.Ramp ramFlo(
        height=0.75*m_flow_nominal,
        duration=3600,
        offset=0.25*m_flow_nominal,
        startTime=120)
        annotation (Placement(transformation(extent={{-90,20},{-70,40}})));
      Modelica.Blocks.Sources.BooleanConstant onSig
        annotation (Placement(transformation(extent={{-66,-50},{-46,-30}})));
      Buildings.Fluid.MixingVolumes.MixingVolume vol(
        m_flow_nominal=m_flow_nominal,
        V=500,                                       nPorts=1, redeclare
          package
          Medium =                                                                        Medium)
        annotation (Placement(transformation(extent={{40,10},{60,30}})));

    equation
      connect(flo.ports[1], ducGUV.port_a)
        annotation (Line(points={{-28,0},{-10,0}}, color={0,127,255}));
      connect(onSig.y, ducGUV.u) annotation (Line(points={{-45,-40},{-20,-40},{-20,
              -8},{-12,-8}}, color={255,0,255}));
      connect(ramFlo.y, flo.m_flow_in) annotation (Line(points={{-69,30},{-58,30},{-58,
              8},{-50,8}}, color={0,0,127}));
      connect(ducGUV.port_b, vol.ports[1]) annotation (Line(points={{10,0},{50,0},
              {50,10}},                    color={0,127,255}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=3840, __Dymola_Algorithm="Dassl"));
    end DuctGUVExample;

    model GUVExample
      replaceable package Medium = Buildings.Media.Air(extraPropertiesNames={"CO2", "SARS-CoV-2"}) "Medium model for air";
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal=1.2*vol.V/3600
        "Design mass flow rate";
      Buildings.Fluid.MixingVolumes.MixingVolumeMoistAir
                                               vol(
        redeclare package Medium = Medium,
        final energyDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial,
        final massDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial,
        final V=500,
        final C_start={0,0},
        final mSenFac=1,
        final m_flow_nominal=m_flow_nominal,
        final prescribedHeatFlowRate=true,
        final nPorts=1,
        m_flow_small=1E-4*abs(m_flow_nominal),
        allowFlowReversal=true,
        final use_C_flow=true)        "Room air volume"
        annotation (Placement(transformation(extent={{38,-10},{58,10}})));
      Modelica.Blocks.Sources.Constant const(k=0)
        annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
      Modelica.Blocks.Sources.RealExpression speConc2(y=vol.C[2])
        "SARS-CoV-2 concentration"
        annotation (Placement(transformation(extent={{-90,-50},{-70,-30}})));
      Buildings.Fluid.Sources.MassFlowSource_T sup(
        use_C_in=false,
        C={0.1,0.1},
        use_m_flow_in=false,
        redeclare package Medium = Medium,
        m_flow=m_flow_nominal,
        nPorts=1)
        annotation (Placement(transformation(extent={{-68,-10},{-48,10}})));
      Modelica.Blocks.Sources.BooleanConstant booleanConstant
        annotation (Placement(transformation(extent={{-90,-72},{-70,-52}})));
      Modelica.Blocks.Sources.RealExpression speConc1(y=vol.C[1])
        "CO2 concentration"
        annotation (Placement(transformation(extent={{-90,-32},{-70,-12}})));
      Buildings.AirCleaning.RoomGUV gUV(redeclare package Medium = Medium,
        frad=0.2,
        Eavg=0.5,
        krad={0,0.28},
        kpow=60,
        V=vol.V)
        annotation (Placement(transformation(extent={{-30,-56},{-10,-36}})));
    equation
      connect(const.y, vol.mWat_flow) annotation (Line(points={{-19,50},{28,50},{28,
              8},{36,8}},  color={0,0,127}));
      connect(sup.ports[1], vol.ports[1]) annotation (Line(points={{-48,0},{-4,0},{-4,
              -12},{34,-12},{34,-14},{48,-14},{48,-10}},
                                         color={0,127,255}));
      connect(booleanConstant.y, gUV.u) annotation (Line(points={{-69,-62},{-38,-62},
              {-38,-48},{-32,-48}}, color={255,0,255}));
      connect(speConc1.y, gUV.C[1]) annotation (Line(points={{-69,-22},{-46,-22},{-46,
              -42},{-32,-42}},                         color={0,0,127}));
      connect(speConc2.y, gUV.C[2]) annotation (Line(points={{-69,-40},{-46,-40},{-46,
              -42},{-32,-42}}, color={0,0,127}));
      connect(gUV.yC_flow[1], vol.C_flow[1]) annotation (Line(points={{-9,-42},{0,-42},
              {0,-6},{36,-6}},     color={0,0,127}));
      connect(gUV.yC_flow[2], vol.C_flow[2]) annotation (Line(points={{-9,-42},{8,-42},
              {8,-6},{36,-6}},     color={0,0,127}));
      connect(gUV.diss, vol.heatPort) annotation (Line(points={{-9.6,-54},{4,-54},{4,
              0},{38,0}},    color={191,0,0}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=3600, __Dymola_Algorithm="Dassl"),
        Documentation(info="<html>
<p>This example simulates <a href=\"modelica://Buildings.AirCleaning.PAC\">AirCleaning.GUV</a> for a scenario with CO<sub>2. </sub>and SARS-CoV-2 as trace species. Airflow with a constant rate and concentration of trace species is supplied to a mixing volume. A GUV device is connected to the mixing volume to simulate removal of trace species from the volume, where the GUV does not inactivate CO<sub>2</sub> but inactivates SARS-CoV-2.</p>
</html>"));
    end GUVExample;
  end Examples;

  annotation ();
end AirCleaning;
