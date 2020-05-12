within Buildings.ThermalZones.EnergyPlus.Validation;
model OneZoneCompareZoneTemperatureOutput
  "This example tests whether the zone mean air temperature is reported correctly as an EnergyPlus output"
  extends Buildings.ThermalZones.EnergyPlus.Validation.OneZone;

  Buildings.ThermalZones.EnergyPlus.OutputVariable zonMeaAirTem(
    key="LIVING ZONE",
    name="Zone Mean Air Temperature",
    y(final unit="K",
      displayUnit="degC"))
    "Block that reads output from EnergyPlus"
    annotation (Placement(transformation(extent={{70,50},{90,70}})));

  annotation (Documentation(info="<html>
<p>
Simple test case that verifies whether the zone mean air temperature is reported correctly by EnergyPlus.
Note that Modelica solves the differential equation for this variable, but this test case
obtains its value from EnergyPlus.
</p>
</html>", revisions="<html>
<ul><li>
April 2, 2020, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"),
 __Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/ThermalZones/EnergyPlus/Validation/OneZoneCompareZoneTemperatureOutput.mos"
        "Simulate and plot"),
experiment(
      StopTime=172800,
      Tolerance=1e-06));
end OneZoneCompareZoneTemperatureOutput;
