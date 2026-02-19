within Buildings.AirCleaning;
model PAC_new
  extends Buildings.Fluid.Interfaces.PartialTwoPortInterface(m_flow_small=0,
    allowFlowReversal = false);

  parameter Modelica.Units.SI.MassFlowRate m_flow_in(min=0)
    "Prescribed mass flow rate through PAC";
  parameter Real eta(min=0, max=1) = 0.9
    "Contaminant removal efficiency (0 = no removal, 1 = full removal)";

  replaceable package Medium = Buildings.Media.Air;

//protected
  //Medium.ThermodynamicState sta;
  //Modelica.Units.SI.MassFraction C_in;
  //Modelica.Units.SI.MassFraction C_out;

equation
  m_flow = m_flow_in;
  //dp = 0;
  // Pressure drop in design flow direction
  dp = port_a.p - port_b.p;

  // Design direction of mass flow rate
  m_flow = port_a.m_flow;
  assert(m_flow > -m_flow_small or allowFlowReversal,
      "Reverting flow occurs even though allowFlowReversal is false");

  // Mass balance (no storage)
  port_a.m_flow + port_b.m_flow = 0;

  // mass balance
  //port_b.m_flow = m_flow;
  // momentum equation (no pressure loss)
  //port_a.p = port_b.p;
  // isenthalpic state transformation (no storage and no loss of energy)
  port_a.h_outflow = if allowFlowReversal then inStream(port_b.h_outflow) else Medium.h_default;
  port_b.h_outflow = inStream(port_a.h_outflow);
  port_a.Xi_outflow = if allowFlowReversal then inStream(port_b.Xi_outflow) else Medium.X_default[1:Medium.nXi];
  port_b.Xi_outflow = inStream(port_a.Xi_outflow);
  port_a.C_outflow = if allowFlowReversal then inStream(port_b.C_outflow) else zeros(Medium.nC);
  port_b.C_outflow = inStream(port_a.C_outflow);
end PAC_new;
