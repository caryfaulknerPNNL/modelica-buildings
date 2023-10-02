within Buildings.Experimental.DHC.Loads.HotWater.Examples;
model HeaterAndFixtureDirectHeatExchangerWithElectricHeat
  extends Modelica.Icons.Example;
  extends BaseClasses.PartialHeaterAndFixture(
    souCol(nPorts=2),
    souDis(nPorts=1),
    sinDis(nPorts=1));
  DirectHeatExchangerWithElectricHeat gen(
    redeclare package Medium = Medium,
    mHotSou_flow_nominal=mHotSou_flow_nominal,
      mDis_flow_nominal=mDis_flow_nominal)
    "Hot water generated by heat exchanger and auxiliary heat"
    annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
equation

  connect(gen.port_a1, souCol.ports[2]) annotation (Line(points={{-50,6},{-54,6},
          {-54,-20},{10,-20},{10,-40}}, color={0,127,255}));
  connect(gen.port_a2, souDis.ports[1]) annotation (Line(points={{-30,-6},{-20,
          -6},{-20,-32},{-30,-32},{-30,-40}}, color={0,127,255}));
  connect(gen.port_b2, sinDis.ports[1])
    annotation (Line(points={{-50,-6},{-60,-6},{-60,-40}}, color={0,127,255}));
  connect(conTSetHotSou.y, gen.TSetHotSou)
    annotation (Line(points={{-69,0},{-51,0}}, color={0,0,127}));
  connect(gen.PHea, PEle) annotation (Line(points={{-29,0},{-20,0},{-20,80},{
          110,80}}, color={0,0,127}));
  connect(tmv.port_hotSou, gen.port_b1) annotation (Line(points={{0,4},{-14,4},
          {-14,6},{-30,6}}, color={0,127,255}));
  annotation (experiment(
      StopTime=86400,
      Tolerance=1e-06),
      __Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/Experimental/DHC/Loads/HotWater/Examples/HeaterAndFixtureDirectHeatExchangerWithElectricHeat.mos"
        "Simulate and plot"),Documentation(info="<html>
<p>
This model implements an example hot water system where the hot water is 
produced using 
The hydronic arrangement modeled in
<a href=\"modelica://Buildings.Experimental.DHC.Loads.HotWater.DirectHeatExchangerWithElectricHeat\">
Buildings.Experimental.DHC.Loads.HotWater.DirectHeatExchangerWithElectricHeat</a>, 
thermostatically mixed down to a distribution temperature, and supplied to a fixture load
defined by a schedule.
</p>
<p>
Such distribution is based on the
<i>Advanced Energy Design Guide for Multifamily Buildings-Achieving Zero Energy</i>
published by ASHRAE in 2022 at <a href=\"https://www.ashrae.org/technical-resources/aedgs/zero-energy-aedg-free-download\">
https://www.ashrae.org/technical-resources/aedgs/zero-energy-aedg-free-download</a>.
</p>
<p align=\"center\">
<img alt=\"image\" src=\"modelica://Buildings/Resources/Images/Experimental/DHC/Loads/HotWater/Example_DirectHeatExchangerWithElectricHeat.png\"/>
</p>
</html>", revisions="<html>
<ul>
<li>
September 11, 2023 by David Blum:<br/>
Extended base class and updated for release.
</li>
<li>
October 20, 2022 by Dre Helmns:<br/>
Initial implementation.
</li>
</ul>
</html>"));
end HeaterAndFixtureDirectHeatExchangerWithElectricHeat;
