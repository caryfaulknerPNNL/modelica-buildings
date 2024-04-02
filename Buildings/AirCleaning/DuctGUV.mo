within Buildings.AirCleaning;
model DuctGUV "In Duct GUV"
  extends Buildings.AirCleaning.BaseClasses.PartialDuctGUV(final
      m_flow_turbulent=if computeFlowResistance then deltaM*m_flow_nominal_pos
         else 0, vol(nPorts=2));

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
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prePow(final alpha=0) if
    addPowerToMedium
    "Prescribed power (=heat and flow work) flow for dynamic model"
    annotation (Placement(transformation(extent={{-20,-60},{0,-40}})));
initial equation
 if computeFlowResistance then
   assert(m_flow_turbulent > 0, "m_flow_turbulent must be bigger than zero.");
 end if;

 assert(m_flow_nominal_pos > 0, "m_flow_nominal_pos must be non-zero. Check parameters.");
equation

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
<p>This model is designed to simulate the inactivation of trace species by in-Duct germicidal ultraviolet (DuctGUV) devices when the DuctGUV is enabled <span style=\"font-family: Courier New;\">uDucGUVEna</span>. The DuctGUV model calculates the concentration of airborne pathogens in the supply air  after using the in-duct GUV device based on the equation described below. The DuctGUV device consumes energy based on a constant power rating, <span style=\"font-family: Courier New;\">kpow</span> whcih is dissipated into the airflow as heat. </p>
<h4>Main Equations </h4>
<p>c<sub>out</sub> = [1-eff<sub>GUV</sub>(flow<sub>GUV</sub>)]c<sub>in</sub> </p>
<p>where c<sub>out</sub> is the pathogen concentration exiting the flow past the device, eff<sub>GUV</sub> is the removal efficency of the pathogen by the PAC, flow<sub>GUV</sub> is the volumetric airflow rate of air through th GUV, c<sub>in</sub> is the pathogen concentration of the flow entering the irradiated space of the device.eff<sub>GUV</sub>*(flow<sub>GUV</sub>) is the pathogen inactivation of the GUV device as a function of the flow rate.

The pathogen inactivation efficency as a function of flow rate is modeled based on the dose of GUV outlined in the following equations

<p>eff<sub>GUV</sub>(flow<sub>GUV</sub>) = 1 -exp(-E<sub>avg</sub>k<sub>rad</sub>V<sub>rad</sub>/flow<sub>GUV</sub>)  </p>

where V<sub>rad</sub> is a volume of irradiated space in the duct, E<sub>avg</sub> is the average fluence rate of the GUV device, and k<sub>rad</sub> is susceptibility of the trace species to the GUV irrdiation,
<h4>Assumptions</h4>
<p>This assumes that the zone is well mixed and the parameters E<sub>avg</sub> and k<sub>rad</sub> are constant for the given simulation when the GUV device is on. </p>
<h4>Important parameters</h4>
<h4>Notes</h4>
<h4>Implementation</h4>
</html>", revisions="<html>
<ul>
<li>
March 29, 2024 by Cary Faulkner:<br/>
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
