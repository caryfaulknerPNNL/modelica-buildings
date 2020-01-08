within Buildings.Fluid.HeatExchangers.CoolingTowers.BaseClasses.Functions.Validation;
model EquivalentHeatCapacity
  "Validation model for the equivalent heat capacity calculation"
  extends Modelica.Icons.Example;

  Real deltaT
    "Change in temperature of the air entering and leaving the cooling tower";

  Modelica.SIunits.Temperature TIn[1,5]=
    [283.15, 288.15, 293.15, 298.15, 303.15]
     "Inlet temperatures";

  Modelica.SIunits.Temperature TOut[1,5]=
    [TIn[1,1]+deltaT,TIn[1,2]+deltaT, TIn[1,3]+deltaT, TIn[1,4]+deltaT, TIn[1,5]+deltaT]
     "Outlet temperatures";

  Modelica.SIunits.SpecificHeatCapacity cpe10
    "Equivalent specific heat capacity with 10degC inlet temperature";
  Modelica.SIunits.SpecificHeatCapacity cpe15
    "Equivalent specific heat capacity with 15degC inlet temperature";
  Modelica.SIunits.SpecificHeatCapacity cpe20
    "Equivalent specific heat capacity with 20degC inlet temperature";
  Modelica.SIunits.SpecificHeatCapacity cpe25
    "Equivalent specific heat capacity with 25degC inlet temperature";
  Modelica.SIunits.SpecificHeatCapacity cpe30
    "Equivalent specific heat capacity with 30degC inlet temperature";

equation
  deltaT = time;
  cpe10 = Buildings.Fluid.HeatExchangers.CoolingTowers.BaseClasses.Functions.equivalentHeatCapacity(
    TIn = TIn[1,1], TOut = TOut[1,1]);
  cpe15 = Buildings.Fluid.HeatExchangers.CoolingTowers.BaseClasses.Functions.equivalentHeatCapacity(
    TIn = TIn[1,2], TOut = TOut[1,2]);
  cpe20 = Buildings.Fluid.HeatExchangers.CoolingTowers.BaseClasses.Functions.equivalentHeatCapacity(
    TIn = TIn[1,3], TOut = TOut[1,3]);
  cpe25 = Buildings.Fluid.HeatExchangers.CoolingTowers.BaseClasses.Functions.equivalentHeatCapacity(
    TIn = TIn[1,4], TOut = TOut[1,4]);
  cpe30 = Buildings.Fluid.HeatExchangers.CoolingTowers.BaseClasses.Functions.equivalentHeatCapacity(
    TIn = TIn[1,5], TOut = TOut[1,5]);

  annotation (
    experiment(StartTime=0, Tolerance=1e-06, StopTime=50),
    __Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/Fluid/HeatExchangers/CoolingTowers/BaseClasses/Functions/Validation/EquivalentHeatCapacity.mos"
        "Simulate and plot"),
  Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>This model validates the computation of the equivalent heat capacity for five inlet temperature conditions and variable changes in temperature between inlet and outlet airflows. </p>
</html>", revisions="<html>
<ul>
<li>
January 6, 2020, by Kathryn Hinkelman:<br/>
First implementation.<br/>
</li>
</ul>
</html>"));
end EquivalentHeatCapacity;
