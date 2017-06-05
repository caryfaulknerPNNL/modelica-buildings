within Buildings.Experimental.OpenBuildingControl.CDL.Continuous;
block SignalRanker "Ranks output signals such that y[i] >= y[i+1]"
  parameter Integer nin=1 "Number of inputs";
  final parameter Integer nout=nin "Number of outputs";
  Interfaces.RealInput u[nin] "Connector of Real input signals"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
  Interfaces.RealOutput y[nout] "Connector of Real output signals"
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));

protected
  Real t "Temporary variable";
algorithm
  y[:] := u[:];
  for i in 1:nin loop
    for j in 1:nin-1 loop
    if y[j] < y[j+1] then
      t      := y[j+1];
      y[j+1] := y[j];
      y[j]   := t;
    end if;
   end for;
  end for;
  annotation (
defaultComponentName="sigRan",
Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
          100,100}}), graphics={Rectangle(
        extent={{-100,-100},{100,100}},
        lineColor={0,0,127},
        fillColor={255,255,255},
        fillPattern=FillPattern.Solid), Text(
        extent={{-150,150},{150,110}},
        textString="%name",
        lineColor={0,0,255}),
               Text(
          extent={{-94,34},{96,-164}},
          lineColor={0,0,255},
          textString="y[i] >= y[i+1]")}),
Documentation(info="<html>
<p>
Block that sorts the input signal <code>u[:]</code> such that the output
signal satisfies <code>y[i] &gt;= y[i+1]</code> for all <code>i=1, ..., nin-1</code>.
</p>
<p>
This block may for example be used in a variable air volume flow
controller to access the position of the dampers that are most open.
</p>
</html>",
revisions="<html>
<ul>
<li>
January 10, 2017, by Milica Grahovac:<br/>
Initial CDL implementation.
</li>
<li>
November 21, 2011, by Michael Wetter:<br/>
Removed <code>assert</code> statement.
</li>
<li>
November 25, 2008, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"));
end SignalRanker;
