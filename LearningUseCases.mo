package LearningUseCases
  package E1
    model simplepipepump
      Modelica.Fluid.Pipes.StaticPipe pipe(length = 0.5, diameter = 0.254, redeclare package Medium = Modelica.Media.Examples.TwoPhaseWater, nParallel = 1, height_ab = 0.5, FlowModel = 0.1)  annotation(
        Placement(transformation(origin = {1, 38}, extent = {{-35, -14}, {35, 14}})));
  Modelica.Fluid.Valves.ValveIncompressible valveIncompressible(redeclare package Medium = Modelica.Media.Examples.TwoPhaseWater, dp_nominal = 1e4, m_flow_nominal = 1, filteredOpening = false, redeclare function valveCharacteristic = Modelica.Fluid.Valves.BaseClasses.ValveCharacteristics.linear)  annotation(
        Placement(transformation(origin = {60, 38}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Machines.Pump pump(redeclare package Medium = Modelica.Media.Examples.TwoPhaseWater, N_nominal = 100, nParallel = 1, use_powerCharacteristic = true, redeclare function powerCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticPower, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.linearFlow, V = 2)  annotation(
        Placement(transformation(origin = {-60, 38}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Modelica.Media.Examples.TwoPhaseWater, V = 2, nPorts = 1)  annotation(
        Placement(transformation(origin = {-100, 36}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(pipe.port_b, valveIncompressible.port_a) annotation(
        Line(points = {{36, 38}, {50, 38}}, color = {0, 127, 255}));
      connect(pipe.port_a, pump.port_b) annotation(
        Line(points = {{-34, 38}, {-50, 38}}, color = {0, 127, 255}));
  connect(volume.ports[1], pump.port_a) annotation(
        Line(points = {{-100, 26}, {-70, 26}, {-70, 38}}, color = {0, 127, 255}));
    end simplepipepump;
  end E1;
  annotation(
    uses(Modelica(version = "4.0.0")));
end LearningUseCases;
