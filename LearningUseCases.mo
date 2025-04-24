package LearningUseCases
  package E1
    model simplepipepump
      Modelica.Fluid.Pipes.StaticPipe pipe(length = 0.5, diameter = 0.254, redeclare package Medium = Modelica.Media.Examples.TwoPhaseWater, nParallel = 1, height_ab = 0.5, FlowModel = 0.1) annotation(
        Placement(transformation(origin = {1, 38}, extent = {{-35, -14}, {35, 14}})));
      Modelica.Fluid.Valves.ValveIncompressible valveIncompressible(redeclare package Medium = Modelica.Media.Examples.TwoPhaseWater, dp_nominal = 1e4, m_flow_nominal = 1, filteredOpening = false, redeclare function valveCharacteristic = Modelica.Fluid.Valves.BaseClasses.ValveCharacteristics.linear) annotation(
        Placement(transformation(origin = {60, 38}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Fluid.Machines.Pump pump(redeclare package Medium = Modelica.Media.Examples.TwoPhaseWater, N_nominal = 100, nParallel = 1, use_powerCharacteristic = true, redeclare function powerCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticPower, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.linearFlow, V = 2) annotation(
        Placement(transformation(origin = {-60, 38}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Modelica.Media.Examples.TwoPhaseWater, V = 2, nPorts = 1) annotation(
        Placement(transformation(origin = {-100, 36}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(pipe.port_b, valveIncompressible.port_a) annotation(
        Line(points = {{36, 38}, {50, 38}}, color = {0, 127, 255}));
      connect(pipe.port_a, pump.port_b) annotation(
        Line(points = {{-34, 38}, {-50, 38}}, color = {0, 127, 255}));
      connect(volume.ports[1], pump.port_a) annotation(
        Line(points = {{-100, 26}, {-70, 26}, {-70, 38}}, color = {0, 127, 255}));
    end simplepipepump;

    model pumpingOfFluidFromPointAtoB
      Modelica.Fluid.Vessels.OpenTank tank(nPorts = 1, height = 1, crossArea = 0.1, redeclare package Medium = Modelica.Media.Examples.TwoPhaseWater, use_portsData = false, final portsData) annotation(
        Placement(transformation(origin = {-128, -2}, extent = {{-20, -20}, {20, 20}})));
      Modelica.Fluid.Vessels.OpenTank tank1(nPorts = 1, height = 1, crossArea = 0.1, redeclare package Medium = Modelica.Media.Examples.TwoPhaseWater) annotation(
        Placement(transformation(origin = {74, 70}, extent = {{-20, -20}, {20, 20}})));
      Modelica.Fluid.Machines.ControlledPump pump(redeclare package Medium = Modelica.Media.Examples.TwoPhaseWater, nParallel = 1, p_a_nominal = 1e5, p_b_nominal = 5e4, m_flow_nominal = 0.1, control_m_flow = true, V = 0.1) annotation(
        Placement(transformation(origin = {-58, 0}, extent = {{-18, -18}, {18, 18}})));
    equation
      connect(tank.ports[1], pump.port_a) annotation(
        Line(points = {{-128, -22}, {-127, -22}, {-127, -42}, {-130, -42}, {-130, -64}, {-76, -64}, {-76, 0}}, color = {0, 127, 255}));
      connect(pump.port_b, tank1.ports[1]) annotation(
        Line(points = {{-40, 0}, {74, 0}, {74, 50}}, color = {0, 127, 255}));
    end pumpingOfFluidFromPointAtoB;

    model PumpingSystem "Model of a pumping system for drinking water"
      extends Modelica.Icons.Example;
      replaceable package Medium = Modelica.Media.Water.StandardWaterOnePhase constrainedby Modelica.Media.Interfaces.PartialMedium;
      //replaceable package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater
      //  constrainedby Modelica.Media.Interfaces.PartialMedium;
      Modelica.Fluid.Sources.FixedBoundary source(nPorts = 1, use_T = true, T = Modelica.Units.Conversions.from_degC(20), p = system.p_ambient, redeclare package Medium = Modelica.Media.Examples.TwoPhaseWater) annotation(
        Placement(transformation(extent = {{-100, -80}, {-80, -60}})));
      Modelica.Fluid.Pipes.StaticPipe pipe(allowFlowReversal = true, length = 100, height_ab = 50, diameter = 0.3, redeclare package Medium = Modelica.Media.Examples.TwoPhaseWater) annotation(
        Placement(transformation(origin = {-20, -51}, extent = {{-9, -10}, {11, 10}}, rotation = 90)));
      Modelica.Fluid.Machines.PrescribedPump pumps(checkValve = true, checkValveHomotopy = Modelica.Fluid.Types.CheckValveHomotopyType.Closed, N_nominal = 1200, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), use_N_in = true, nParallel = 1, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, V(displayUnit = "l") = 0.05, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, redeclare package Medium = Modelica.Media.Examples.TwoPhaseWater, p_b_start = 600000, T_start = system.T_start) annotation(
        Placement(transformation(extent = {{-68, -80}, {-48, -60}})));
      Modelica.Fluid.Valves.ValveLinear userValve(allowFlowReversal = false, dp_nominal = 200000, m_flow_nominal = 400, redeclare package Medium = Modelica.Media.Examples.TwoPhaseWater) annotation(
        Placement(transformation(extent = {{58, -38}, {74, -22}})));
      Modelica.Fluid.Sources.FixedBoundary sink(p = system.p_ambient, T = system.T_ambient, nPorts = 2, redeclare package Medium = Modelica.Media.Examples.TwoPhaseWater) annotation(
        Placement(transformation(extent = {{100, -40}, {80, -20}})));
      inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial) annotation(
        Placement(transformation(extent = {{60, -96}, {80, -76}})));
    
      model pumpSpeedInput
      
      parameter Real n_speed = 1000 ;
    Modelica.Blocks.Interfaces.RealOutput y annotation(
          Placement(transformation(origin = {-144, 0}, extent = {{100, -10}, {120, 10}}), iconTransformation(origin = {-98, -10}, extent = {{100, -10}, {120, 10}})));
      equation
    n_speed =  y ;
      annotation(
          Icon(graphics = {Rectangle(origin = {-23, -6}, lineColor = {255, 170, 127}, fillColor = {85, 255, 127}, fillPattern = FillPattern.Solid, extent = {{-25, 42}, {25, -42}})}));
    end pumpSpeedInput;
    
      pumpSpeedInput pumpSpeedInput1 annotation(
        Placement(transformation(origin = {-56, -46}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    pumpSpeedInput valueOpeningConstant(n_speed = 0.5)  annotation(
        Placement(transformation(origin = {68, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    equation
      connect(userValve.port_b, sink.ports[1]) annotation(
        Line(points = {{74, -30}, {77, -30}, {77, -28}, {80, -28}}, color = {0, 127, 255}));
      connect(source.ports[1], pumps.port_a) annotation(
        Line(points = {{-80, -70}, {-74, -70}, {-68, -70}}, color = {0, 127, 255}));
      connect(pipe.port_a, pumps.port_b) annotation(
        Line(points = {{-20, -60}, {-20, -70}, {-48, -70}}, color = {0, 127, 255}));
      connect(pipe.port_b, userValve.port_a) annotation(
        Line(points = {{-20, -40}, {13, -40}, {13, -36}, {12, -36}, {12, -30}, {58, -30}}, color = {0, 127, 255}));
      connect(pumpSpeedInput1.y, pumps.N_in) annotation(
        Line(points = {{-57, -47}, {-57, -60}, {-58, -60}}, color = {0, 0, 127}));
    connect(valueOpeningConstant.y, userValve.opening) annotation(
        Line(points = {{67, -7}, {67, -12}, {66, -12}, {66, -24}}, color = {0, 0, 127}));
      annotation(
        Documentation(info = "<html>
    <p>
    Water is pumped from a source by a pump (fitted with check valves), through a pipe whose outlet is 50 m higher than the source, into a reservoir. The users are represented by an equivalent valve, connected to the reservoir.
    </p>
    <p>
    The water controller is a simple on-off controller, regulating on the gauge pressure measured at the base of the tower; the output of the controller is the rotational speed of the pump, which is represented by the output of a first-order system. A small but nonzero rotational speed is used to represent the standby state of the pumps, in order to avoid singularities in the flow characteristic.
    </p>
    <p>When the simulation starts, the level is above the set point, so the initial state of the pump controller is off. Hence, the check valve of the pump is engaged. In order to facilitate the solution of the initialization problem, the <code>homotopyType</code> parameter is set accordingly.
    </p>
    <p>
    Simulate for 2000 s. When the valve is opened at time t=200, the pump starts turning on and off to keep the reservoir level around 2 meters, which roughly corresponds to a gauge pressure of 200 mbar.
    </p>
    
    <img src=\"modelica://Modelica/Resources/Images/Fluid/Examples/PumpingSystem.png\" border=\"1\"
         alt=\"PumpingSystem.png\">
    </html>", revisions = "<html>
    <ul>
    <li><em>Jan 2009</em>
        by R&uuml;diger Franke:<br>
           Reduce diameters of pipe and reservoir ports; use separate port for measurement of reservoirPressure, avoiding disturbances due to pressure losses.</li>
    <li><em>1 Oct 2007</em>
        by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
           Parameters updated.</li>
    <li><em>2 Nov 2005</em>
        by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
           Created.</li>
    </ul>
    </html>"),
        experiment(StopTime = 2000, Interval = 0.4, Tolerance = 1e-006),
        Diagram);
    end PumpingSystem;
  end E1;
  annotation(
    uses(Modelica(version = "4.0.0")));
end LearningUseCases;
