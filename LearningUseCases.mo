package LearningUseCases
  package E1
    model simplepipepump
  Modelica.Fluid.Sources.FixedBoundary source(redeclare package Medium = Medium, T = Modelica.Units.Conversions.from_degC(20), nPorts = 1, p = system.p_ambient, use_T = true) annotation(
        Placement(transformation(extent = {{-100, -80}, {-80, -60}})));
  Modelica.Fluid.Pipes.StaticPipe pipe(redeclare package Medium = Medium, allowFlowReversal = true, diameter = 0.3, height_ab = 50, length = 100) annotation(
        Placement(transformation(origin = {-30, -51}, extent = {{-9, -10}, {11, 10}}, rotation = 90)));
  Modelica.Fluid.Machines.PrescribedPump pumps(redeclare package Medium = Medium, N_nominal = 1200, T_start = system.T_start, V(displayUnit = "l") = 0.05, checkValve = true, checkValveHomotopy = Modelica.Fluid.Types.CheckValveHomotopyType.Closed, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nParallel = 1, p_b_start = 600000, use_N_in = true) annotation(
        Placement(transformation(extent = {{-68, -80}, {-48, -60}})));
  Modelica.Fluid.Vessels.OpenTank reservoir(redeclare package Medium = Medium, T_start = Modelica.Units.Conversions.from_degC(20), crossArea = 50, height = 3, level_start = 2.2, nPorts = 3, portsData = {Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.3), Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.3), Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.01)}, use_portsData = true) annotation(
        Placement(transformation(extent = {{-20, -16}, {0, 4}})));
  Modelica.Fluid.Valves.ValveLinear userValve(redeclare package Medium = Medium, allowFlowReversal = false, dp_nominal = 200000, m_flow_nominal = 400) annotation(
        Placement(transformation(extent = {{58, -38}, {74, -22}})));
  Modelica.Fluid.Sources.FixedBoundary sink(redeclare package Medium = Medium, T = system.T_ambient, nPorts = 2, p = system.p_ambient) annotation(
        Placement(transformation(extent = {{100, -40}, {80, -20}})));
  Modelica.Blocks.Sources.Step valveOpening(offset = 1e-6, startTime = 200) annotation(
        Placement(transformation(extent = {{56, 0}, {76, 20}})));
  Modelica.Blocks.Sources.Constant RelativePressureSetPoint(k = 2e4) annotation(
        Placement(transformation(extent = {{-100, 60}, {-80, 80}})));
  Modelica.Blocks.Logical.OnOffController controller(bandwidth = 4000, pre_y_start = false) annotation(
        Placement(transformation(extent = {{-40, 60}, {-20, 80}})));
  Modelica.Blocks.Logical.TriggeredTrapezoid PumpRPMGenerator(amplitude = 1200, falling = 3, offset = 0.001, rising = 3) annotation(
        Placement(transformation(extent = {{0, 60}, {20, 80}})));
  Modelica.Fluid.Sensors.RelativePressure reservoirPressure(redeclare package Medium = Medium) annotation(
        Placement(transformation(extent = {{10, -12}, {30, -32}})));
  Modelica.Blocks.Continuous.FirstOrder PT1(T = 2, initType = Modelica.Blocks.Types.Init.InitialState, y_start = 0) annotation(
        Placement(transformation(extent = {{40, 60}, {60, 80}})));
  Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial) annotation(
        Placement(transformation(extent = {{60, -96}, {80, -76}})));
  Modelica.Fluid.Sources.FixedBoundary source1(redeclare package Medium = Medium, T = Modelica.Units.Conversions.from_degC(20), nPorts = 1, p = system1.p_ambient, use_T = true) annotation(
        Placement(transformation(extent = {{-100, -80}, {-80, -60}})));
  Modelica.Fluid.Pipes.StaticPipe pipe1(redeclare package Medium = Medium, allowFlowReversal = true, diameter = 0.3, height_ab = 50, length = 100) annotation(
        Placement(transformation(origin = {-30, -51}, extent = {{-9, -10}, {11, 10}}, rotation = 90)));
  Modelica.Fluid.Machines.PrescribedPump pumps1(redeclare package Medium = Medium, N_nominal = 1200, T_start = system1.T_start, V(displayUnit = "l") = 0.05, checkValve = true, checkValveHomotopy = Modelica.Fluid.Types.CheckValveHomotopyType.Closed, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nParallel = 1, p_b_start = 600000, use_N_in = true) annotation(
        Placement(transformation(extent = {{-68, -80}, {-48, -60}})));
  Modelica.Fluid.Vessels.OpenTank reservoir1(redeclare package Medium = Medium, T_start = Modelica.Units.Conversions.from_degC(20), crossArea = 50, height = 3, level_start = 2.2, nPorts = 3, portsData = {Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.3), Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.3), Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.01)}, use_portsData = true) annotation(
        Placement(transformation(extent = {{-20, -16}, {0, 4}})));
  Modelica.Fluid.Valves.ValveLinear userValve1(redeclare package Medium = Medium, allowFlowReversal = false, dp_nominal = 200000, m_flow_nominal = 400) annotation(
        Placement(transformation(extent = {{58, -38}, {74, -22}})));
  Modelica.Fluid.Sources.FixedBoundary sink1(redeclare package Medium = Medium, T = system1.T_ambient, nPorts = 2, p = system1.p_ambient) annotation(
        Placement(transformation(extent = {{100, -40}, {80, -20}})));
  Modelica.Blocks.Sources.Step valveOpening1(offset = 1e-6, startTime = 200) annotation(
        Placement(transformation(extent = {{56, 0}, {76, 20}})));
  Modelica.Blocks.Sources.Constant RelativePressureSetPoint1(k = 2e4) annotation(
        Placement(transformation(extent = {{-100, 60}, {-80, 80}})));
  Modelica.Blocks.Logical.OnOffController controller1(bandwidth = 4000, pre_y_start = false) annotation(
        Placement(transformation(extent = {{-40, 60}, {-20, 80}})));
  Modelica.Blocks.Logical.TriggeredTrapezoid PumpRPMGenerator1(amplitude = 1200, falling = 3, offset = 0.001, rising = 3) annotation(
        Placement(transformation(extent = {{0, 60}, {20, 80}})));
  Modelica.Fluid.Sensors.RelativePressure reservoirPressure1(redeclare package Medium = Medium) annotation(
        Placement(transformation(extent = {{10, -12}, {30, -32}})));
  Modelica.Blocks.Continuous.FirstOrder PT11(T = 2, initType = Modelica.Blocks.Types.Init.InitialState, y_start = 0) annotation(
        Placement(transformation(extent = {{40, 60}, {60, 80}})));
  Modelica.Fluid.System system1(energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial) annotation(
        Placement(transformation(extent = {{60, -96}, {80, -76}})));
    equation
      connect(userValve.port_b, sink.ports[1]) annotation(
        Line(points = {{74, -30}, {77, -30}, {77, -28}, {80, -28}}, color = {0, 127, 255}));
      connect(source.ports[1], pumps.port_a) annotation(
        Line(points = {{-80, -70}, {-74, -70}, {-68, -70}}, color = {0, 127, 255}));
      connect(valveOpening.y, userValve.opening) annotation(
        Line(points = {{77, 10}, {98, 10}, {98, -12}, {66, -12}, {66, -23.6}}, color = {0, 0, 127}));
      connect(RelativePressureSetPoint.y, controller.reference) annotation(
        Line(points = {{-79, 70}, {-60, 70}, {-60, 76}, {-42, 76}}, color = {0, 0, 127}));
      connect(controller.y, PumpRPMGenerator.u) annotation(
        Line(points = {{-19, 70}, {-2, 70}}, color = {255, 0, 255}));
      connect(reservoirPressure.p_rel, controller.u) annotation(
        Line(points = {{20, -13}, {20, 50}, {-52, 50}, {-52, 64}, {-42, 64}}, color = {0, 0, 127}));
      connect(reservoirPressure.port_b, sink.ports[2]) annotation(
        Line(points = {{30, -22}, {44, -22}, {44, -48}, {80, -48}, {80, -32}}, color = {0, 127, 255}, pattern = LinePattern.Dot));
      connect(PumpRPMGenerator.y, PT1.u) annotation(
        Line(points = {{21, 70}, {38, 70}}, color = {0, 0, 127}));
      connect(PT1.y, pumps.N_in) annotation(
        Line(points = {{61, 70}, {74, 70}, {74, 30}, {-58, 30}, {-58, -60}}, color = {0, 0, 127}));
      connect(pipe.port_a, pumps.port_b) annotation(
        Line(points = {{-30, -60}, {-30, -70}, {-48, -70}}, color = {0, 127, 255}));
      connect(reservoir.ports[1], pipe.port_b) annotation(
        Line(points = {{-12.6667, -16}, {-12.6667, -30}, {-30, -30}, {-30, -40}}, color = {0, 127, 255}));
      connect(reservoir.ports[3], reservoirPressure.port_a) annotation(
        Line(points = {{-7.33333, -16}, {-7, -16}, {-7, -22}, {10, -22}}, color = {0, 127, 255}, pattern = LinePattern.Dot));
      connect(reservoir.ports[2], userValve.port_a) annotation(
        Line(points = {{-10, -16}, {-10, -30}, {58, -30}}, color = {0, 127, 255}));
      connect(userValve1.port_b, sink1.ports[1]) annotation(
        Line(points = {{74, -30}, {77, -30}, {77, -28}, {80, -28}}, color = {0, 127, 255}));
      connect(source1.ports[1], pumps1.port_a) annotation(
        Line(points = {{-80, -70}, {-74, -70}, {-68, -70}}, color = {0, 127, 255}));
      connect(valveOpening1.y, userValve1.opening) annotation(
        Line(points = {{77, 10}, {98, 10}, {98, -12}, {66, -12}, {66, -23.6}}, color = {0, 0, 127}));
      connect(RelativePressureSetPoint1.y, controller1.reference) annotation(
        Line(points = {{-79, 70}, {-60, 70}, {-60, 76}, {-42, 76}}, color = {0, 0, 127}));
      connect(controller1.y, PumpRPMGenerator1.u) annotation(
        Line(points = {{-19, 70}, {-2, 70}}, color = {255, 0, 255}));
      connect(reservoirPressure1.p_rel, controller1.u) annotation(
        Line(points = {{20, -13}, {20, 50}, {-52, 50}, {-52, 64}, {-42, 64}}, color = {0, 0, 127}));
      connect(reservoirPressure1.port_b, sink1.ports[2]) annotation(
        Line(points = {{30, -22}, {44, -22}, {44, -48}, {80, -48}, {80, -32}}, color = {0, 127, 255}, pattern = LinePattern.Dot));
      connect(PumpRPMGenerator1.y, PT11.u) annotation(
        Line(points = {{21, 70}, {38, 70}}, color = {0, 0, 127}));
      connect(PT11.y, pumps1.N_in) annotation(
        Line(points = {{61, 70}, {74, 70}, {74, 30}, {-58, 30}, {-58, -60}}, color = {0, 0, 127}));
      connect(pipe1.port_a, pumps1.port_b) annotation(
        Line(points = {{-30, -60}, {-30, -70}, {-48, -70}}, color = {0, 127, 255}));
      connect(reservoir1.ports[1], pipe1.port_b) annotation(
        Line(points = {{-12.6667, -16}, {-12.6667, -30}, {-30, -30}, {-30, -40}}, color = {0, 127, 255}));
      connect(reservoir1.ports[3], reservoirPressure1.port_a) annotation(
        Line(points = {{-7.33333, -16}, {-7, -16}, {-7, -22}, {10, -22}}, color = {0, 127, 255}, pattern = LinePattern.Dot));
      connect(reservoir1.ports[2], userValve1.port_a) annotation(
        Line(points = {{-10, -16}, {-10, -30}, {58, -30}}, color = {0, 127, 255}));
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

    partial package MediumLN2
    
    
      "Base class for two phase medium of one substance"
      extends PartialPureSubstance(redeclare replaceable record FluidConstants =
            Modelica.Media.Interfaces.Types.TwoPhase.FluidConstants);
      constant Boolean smoothModel=false
        "True if the (derived) model should not generate state events";
      constant Boolean onePhase=false
        "True if the (derived) model should never be called with two-phase inputs";
    
      constant FluidConstants[nS] fluidConstants "Constant data for the fluid";
    
      redeclare replaceable record extends ThermodynamicState
        "Thermodynamic state of two phase medium"
        FixedPhase phase(min=0, max=2)
          "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use";
      end ThermodynamicState;
    
      redeclare replaceable partial model extends BaseProperties
        "Base properties (p, d, T, h, u, R_s, MM, sat) of two phase medium"
        SaturationProperties sat "Saturation properties at the medium pressure";
      end BaseProperties;
    
      replaceable partial function setDewState
        "Return the thermodynamic state on the dew line"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation point";
        input FixedPhase phase(
          min=1,
          max=2) = 1 "Phase: default is one phase";
        output ThermodynamicState state "Complete thermodynamic state info";
      end setDewState;
    
      replaceable partial function setBubbleState
        "Return the thermodynamic state on the bubble line"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation point";
        input FixedPhase phase(
          min=1,
          max=2) = 1 "Phase: default is one phase";
        output ThermodynamicState state "Complete thermodynamic state info";
      end setBubbleState;
    
      redeclare replaceable partial function extends setState_dTX
        "Return thermodynamic state as function of d, T and composition X or Xi"
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
      end setState_dTX;
    
      redeclare replaceable partial function extends setState_phX
        "Return thermodynamic state as function of p, h and composition X or Xi"
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
      end setState_phX;
    
      redeclare replaceable partial function extends setState_psX
        "Return thermodynamic state as function of p, s and composition X or Xi"
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
      end setState_psX;
    
      redeclare replaceable partial function extends setState_pTX
        "Return thermodynamic state as function of p, T and composition X or Xi"
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
      end setState_pTX;
    
      replaceable function setSat_T
        "Return saturation property record from temperature"
        extends Modelica.Icons.Function;
        input Temperature T "Temperature";
        output SaturationProperties sat "Saturation property record";
      algorithm
        sat.Tsat := T;
        sat.psat := saturationPressure(T);
      end setSat_T;
    
      replaceable function setSat_p
        "Return saturation property record from pressure"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        output SaturationProperties sat "Saturation property record";
      algorithm
        sat.psat := p;
        sat.Tsat := saturationTemperature(p);
      end setSat_p;
    
      replaceable partial function bubbleEnthalpy
        "Return bubble point specific enthalpy"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output SI.SpecificEnthalpy hl "Boiling curve specific enthalpy";
      end bubbleEnthalpy;
    
      replaceable partial function dewEnthalpy
        "Return dew point specific enthalpy"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output SI.SpecificEnthalpy hv "Dew curve specific enthalpy";
      end dewEnthalpy;
    
      replaceable partial function bubbleEntropy
        "Return bubble point specific entropy"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output SI.SpecificEntropy sl "Boiling curve specific entropy";
      end bubbleEntropy;
    
      replaceable partial function dewEntropy "Return dew point specific entropy"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output SI.SpecificEntropy sv "Dew curve specific entropy";
      end dewEntropy;
    
      replaceable partial function bubbleDensity "Return bubble point density"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output Density dl "Boiling curve density";
      end bubbleDensity;
    
      replaceable partial function dewDensity "Return dew point density"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output Density dv "Dew curve density";
      end dewDensity;
    
      replaceable partial function saturationPressure
        "Return saturation pressure"
        extends Modelica.Icons.Function;
        input Temperature T "Temperature";
        output AbsolutePressure p "Saturation pressure";
      end saturationPressure;
    
      replaceable partial function saturationTemperature
        "Return saturation temperature"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        output Temperature T "Saturation temperature";
      end saturationTemperature;
    
      replaceable function saturationPressure_sat "Return saturation pressure"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output AbsolutePressure p "Saturation pressure";
      algorithm
        p := sat.psat;
      end saturationPressure_sat;
    
      replaceable function saturationTemperature_sat
        "Return saturation temperature"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output Temperature T "Saturation temperature";
      algorithm
        T := sat.Tsat;
      end saturationTemperature_sat;
    
      replaceable partial function saturationTemperature_derp
        "Return derivative of saturation temperature w.r.t. pressure"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        output DerTemperatureByPressure dTp
          "Derivative of saturation temperature w.r.t. pressure";
      end saturationTemperature_derp;
    
      replaceable function saturationTemperature_derp_sat
        "Return derivative of saturation temperature w.r.t. pressure"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output DerTemperatureByPressure dTp
          "Derivative of saturation temperature w.r.t. pressure";
      algorithm
        dTp := saturationTemperature_derp(sat.psat);
      end saturationTemperature_derp_sat;
    
      replaceable partial function surfaceTension
        "Return surface tension sigma in the two phase region"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output SurfaceTension sigma
          "Surface tension sigma in the two phase region";
      end surfaceTension;
    
      redeclare replaceable function extends molarMass
        "Return the molar mass of the medium"
      algorithm
        MM := fluidConstants[1].molarMass;
      end molarMass;
    
      replaceable partial function dBubbleDensity_dPressure
        "Return bubble point density derivative"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output DerDensityByPressure ddldp "Boiling curve density derivative";
      end dBubbleDensity_dPressure;
    
      replaceable partial function dDewDensity_dPressure
        "Return dew point density derivative"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output DerDensityByPressure ddvdp "Saturated steam density derivative";
      end dDewDensity_dPressure;
    
      replaceable partial function dBubbleEnthalpy_dPressure
        "Return bubble point specific enthalpy derivative"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output DerEnthalpyByPressure dhldp
          "Boiling curve specific enthalpy derivative";
      end dBubbleEnthalpy_dPressure;
    
      replaceable partial function dDewEnthalpy_dPressure
        "Return dew point specific enthalpy derivative"
        extends Modelica.Icons.Function;
    
        input SaturationProperties sat "Saturation property record";
        output DerEnthalpyByPressure dhvdp
          "Saturated steam specific enthalpy derivative";
      end dDewEnthalpy_dPressure;
    
      redeclare replaceable function specificEnthalpy_pTX
        "Return specific enthalpy from pressure, temperature and mass fraction"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[:] "Mass fractions";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "Specific enthalpy at p, T, X";
      algorithm
        h := specificEnthalpy(setState_pTX(
                p,
                T,
                X,
                phase));
      end specificEnthalpy_pTX;
    
      redeclare replaceable function temperature_phX
        "Return temperature from p, h, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[:] "Mass fractions";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output Temperature T "Temperature";
      algorithm
        T := temperature(setState_phX(
                p,
                h,
                X,
                phase));
      end temperature_phX;
    
      redeclare replaceable function density_phX
        "Return density from p, h, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[:] "Mass fractions";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        d := density(setState_phX(
                p,
                h,
                X,
                phase));
      end density_phX;
    
      redeclare replaceable function temperature_psX
        "Return temperature from p, s, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:] "Mass fractions";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output Temperature T "Temperature";
      algorithm
        T := temperature(setState_psX(
                p,
                s,
                X,
                phase));
      end temperature_psX;
    
      redeclare replaceable function density_psX
        "Return density from p, s, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:] "Mass fractions";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        d := density(setState_psX(
                p,
                s,
                X,
                phase));
      end density_psX;
    
      redeclare replaceable function specificEnthalpy_psX
        "Return specific enthalpy from p, s, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:] "Mass fractions";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy(setState_psX(
                p,
                s,
                X,
                phase));
      end specificEnthalpy_psX;
    
      redeclare replaceable function setState_pT
        "Return thermodynamic state from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_pTX(
                p,
                T,
                fill(0, 0),
                phase);
      end setState_pT;
    
      redeclare replaceable function setState_ph
        "Return thermodynamic state from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_phX(
                p,
                h,
                fill(0, 0),
                phase);
      end setState_ph;
    
      redeclare replaceable function setState_ps
        "Return thermodynamic state from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_psX(
                p,
                s,
                fill(0, 0),
                phase);
      end setState_ps;
    
      redeclare replaceable function setState_dT
        "Return thermodynamic state from d and T"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_dTX(
                d,
                T,
                fill(0, 0),
                phase);
      end setState_dT;
    
      replaceable function setState_px
        "Return thermodynamic state from pressure and vapour quality"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input MassFraction x "Vapour quality";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_ph(
                p,
                (1 - x)*bubbleEnthalpy(setSat_p(p)) + x*dewEnthalpy(setSat_p(p)),
                2);
      end setState_px;
    
      replaceable function setState_Tx
        "Return thermodynamic state from temperature and vapour quality"
        extends Modelica.Icons.Function;
        input Temperature T "Temperature";
        input MassFraction x "Vapour quality";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_ph(
                saturationPressure_sat(setSat_T(T)),
                (1 - x)*bubbleEnthalpy(setSat_T(T)) + x*dewEnthalpy(setSat_T(T)),
                2);
      end setState_Tx;
    
      replaceable function vapourQuality "Return vapour quality"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output MassFraction x "Vapour quality";
      protected
        constant SpecificEnthalpy eps=1e-8;
      algorithm
        x := min(max((specificEnthalpy(state) - bubbleEnthalpy(setSat_p(pressure(
          state))))/(dewEnthalpy(setSat_p(pressure(state))) - bubbleEnthalpy(
          setSat_p(pressure(state))) + eps), 0), 1);
      end vapourQuality;
    
      redeclare replaceable function density_ph "Return density from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        d := density_phX(
                p,
                h,
                fill(0, 0),
                phase);
      end density_ph;
    
      redeclare replaceable function temperature_ph
        "Return temperature from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output Temperature T "Temperature";
      algorithm
        T := temperature_phX(
                p,
                h,
                fill(0, 0),
                phase);
      end temperature_ph;
    
      redeclare replaceable function pressure_dT "Return pressure from d and T"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output AbsolutePressure p "Pressure";
      algorithm
        p := pressure(setState_dTX(
                d,
                T,
                fill(0, 0),
                phase));
      end pressure_dT;
    
      redeclare replaceable function specificEnthalpy_dT
        "Return specific enthalpy from d and T"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy(setState_dTX(
                d,
                T,
                fill(0, 0),
                phase));
      end specificEnthalpy_dT;
    
      redeclare replaceable function specificEnthalpy_ps
        "Return specific enthalpy from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy_psX(
                p,
                s,
                fill(0, 0));
      end specificEnthalpy_ps;
    
      redeclare replaceable function temperature_ps
        "Return temperature from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output Temperature T "Temperature";
      algorithm
        T := temperature_psX(
                p,
                s,
                fill(0, 0),
                phase);
      end temperature_ps;
    
      redeclare replaceable function density_ps "Return density from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        d := density_psX(
                p,
                s,
                fill(0, 0),
                phase);
      end density_ps;
    
      redeclare replaceable function specificEnthalpy_pT
        "Return specific enthalpy from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy_pTX(
                p,
                T,
                fill(0, 0),
                phase);
      end specificEnthalpy_pT;
    
      redeclare replaceable function density_pT "Return density from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        d := density(setState_pTX(
                p,
                T,
                fill(0, 0),
                phase));
      end density_pT;
    
    
    
    end MediumLN2;
  end E1;
  
  model pumpSpeedInput
  
  parameter Real n_speed = 1000 ;
  Modelica.Blocks.Interfaces.RealOutput y annotation(
      Placement(transformation(origin = {-144, 0}, extent = {{100, -10}, {120, 10}}), iconTransformation(origin = {-98, -10}, extent = {{100, -10}, {120, 10}})));
  equation
  n_speed =  y ;
  annotation(
      Icon(graphics = {Rectangle(origin = {-23, -6}, lineColor = {255, 170, 127}, fillColor = {85, 255, 127}, fillPattern = FillPattern.Solid, extent = {{-25, 42}, {25, -42}})}));
  end pumpSpeedInput;
  
  model PumpingSystem "Model of a pumping system for drinking water"
    extends Modelica.Icons.Example;
    replaceable package Medium = Modelica.Media.Water.StandardWaterOnePhase constrainedby Modelica.Media.Interfaces.PartialMedium;
    //replaceable package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater
    //  constrainedby Modelica.Media.Interfaces.PartialMedium;
    Modelica.Fluid.Sources.FixedBoundary source(nPorts = 1, use_T = true, T = Modelica.Units.Conversions.from_degC(20), p = system.p_ambient, redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{-100, -80}, {-80, -60}})));
    Modelica.Fluid.Pipes.StaticPipe pipe(allowFlowReversal = true, length = 100, height_ab = 50, diameter = 0.3, redeclare package Medium = Medium) annotation(
      Placement(transformation(origin = {-30, -51}, extent = {{-9, -10}, {11, 10}}, rotation = 90)));
    Modelica.Fluid.Machines.PrescribedPump pumps(checkValve = true, checkValveHomotopy = Modelica.Fluid.Types.CheckValveHomotopyType.Closed, N_nominal = 1200, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), use_N_in = true, nParallel = 1, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, V(displayUnit = "l") = 0.05, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, redeclare package Medium = Medium, p_b_start = 600000, T_start = system.T_start) annotation(
      Placement(transformation(extent = {{-68, -80}, {-48, -60}})));
    Modelica.Fluid.Vessels.OpenTank reservoir(T_start = Modelica.Units.Conversions.from_degC(20), use_portsData = true, crossArea = 50, level_start = 2.2, height = 3, nPorts = 3, portsData = {Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.3), Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.3), Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.01)}, redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{-20, -16}, {0, 4}})));
    Modelica.Fluid.Valves.ValveLinear userValve(allowFlowReversal = false, dp_nominal = 200000, m_flow_nominal = 400, redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{58, -38}, {74, -22}})));
    Modelica.Fluid.Sources.FixedBoundary sink(p = system.p_ambient, T = system.T_ambient, nPorts = 2, redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{100, -40}, {80, -20}})));
    inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial) annotation(
      Placement(transformation(extent = {{60, -96}, {80, -76}})));
  pumpSpeedInput pumpSpeedInput1(n_speed = 1200)  annotation(
      Placement(transformation(origin = {-62, -24}, extent = {{-10, -10}, {10, 10}})));
  pumpSpeedInput valveopeningconstant(n_speed = 0.5)  annotation(
      Placement(transformation(origin = {64, 16}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(userValve.port_b, sink.ports[1]) annotation(
      Line(points = {{74, -30}, {80, -30}}, color = {0, 127, 255}));
    connect(source.ports[1], pumps.port_a) annotation(
      Line(points = {{-80, -70}, {-74, -70}, {-68, -70}}, color = {0, 127, 255}));
    connect(pipe.port_a, pumps.port_b) annotation(
      Line(points = {{-30, -60}, {-30, -70}, {-48, -70}}, color = {0, 127, 255}));
    connect(reservoir.ports[1], pipe.port_b) annotation(
      Line(points = {{-12.6667, -16}, {-12.6667, -30}, {-30, -30}, {-30, -40}}, color = {0, 127, 255}));
    connect(reservoir.ports[2], userValve.port_a) annotation(
      Line(points = {{-10, -16}, {-10, -30}, {58, -30}}, color = {0, 127, 255}));
    connect(pumpSpeedInput1.y, pumps.N_in) annotation(
      Line(points = {{-60, -24}, {-58, -24}, {-58, -60}}, color = {0, 0, 127}));
    connect(valveopeningconstant.y, userValve.opening) annotation(
      Line(points = {{66, 16}, {66, -24}}, color = {0, 0, 127}));
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
  
  model PumpingSystem_reservoirbefore "Model of a pumping system for drinking water"
    extends Modelica.Icons.Example;
    replaceable package Medium = Modelica.Media.Water.StandardWaterOnePhase constrainedby Modelica.Media.Interfaces.PartialMedium;
    //replaceable package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater
    //  constrainedby Modelica.Media.Interfaces.PartialMedium;
    Modelica.Fluid.Sources.FixedBoundary source(nPorts = 1, use_T = true, T = Modelica.Units.Conversions.from_degC(20), p = system.p_ambient, redeclare package Medium = Medium) annotation(
      Placement(transformation(origin = {-102, 14}, extent = {{-100, -80}, {-80, -60}})));
    Modelica.Fluid.Pipes.StaticPipe pipe(allowFlowReversal = true, length = 100, height_ab = 50, diameter = 0.3, redeclare package Medium = Medium) annotation(
      Placement(transformation(origin = {4, -51}, extent = {{-9, -10}, {11, 10}})));
    Modelica.Fluid.Machines.PrescribedPump pumps(checkValve = true, checkValveHomotopy = Modelica.Fluid.Types.CheckValveHomotopyType.Closed, N_nominal = 1200, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), use_N_in = true, nParallel = 1, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, V(displayUnit = "l") = 0.05, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, redeclare package Medium = Medium, p_b_start = 600000, T_start = system.T_start) annotation(
      Placement(transformation(extent = {{-68, -80}, {-48, -60}})));
    Modelica.Fluid.Valves.ValveLinear userValve(allowFlowReversal = false, dp_nominal = 200000, m_flow_nominal = 400, redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{58, -38}, {74, -22}})));
    Modelica.Fluid.Sources.FixedBoundary sink(p = system.p_ambient, T = system.T_ambient, nPorts = 2, redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{100, -40}, {80, -20}})));
    inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial) annotation(
      Placement(transformation(extent = {{60, -96}, {80, -76}})));
    pumpSpeedInput pumpSpeedInput1(n_speed = 1200) annotation(
      Placement(transformation(origin = {-56, -44}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    pumpSpeedInput valveopeningconstant(n_speed = 0.5) annotation(
      Placement(transformation(origin = {68, -4}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Vessels.OpenTank reservoir(redeclare package Medium = Medium, T_start = Modelica.Units.Conversions.from_degC(20), crossArea = 50, height = 3, level_start = 2.2, nPorts = 3, portsData = {Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.3), Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.3), Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.01)}, use_portsData = true) annotation(
      Placement(transformation(origin = {-116, -54}, extent = {{-20, -16}, {0, 4}})));
  equation
    connect(userValve.port_b, sink.ports[1]) annotation(
      Line(points = {{74, -30}, {80, -30}}, color = {0, 127, 255}));
    connect(pipe.port_a, pumps.port_b) annotation(
      Line(points = {{-5, -51}, {-5, -70}, {-48, -70}}, color = {0, 127, 255}));
    connect(pumpSpeedInput1.y, pumps.N_in) annotation(
      Line(points = {{-57, -45}, {-57, -60}, {-58, -60}}, color = {0, 0, 127}));
    connect(valveopeningconstant.y, userValve.opening) annotation(
      Line(points = {{69, -5}, {69, -16}, {72, -16}, {72, -15}, {66, -15}, {66, -24}}, color = {0, 0, 127}));
    connect(reservoir.ports[2], pumps.port_a) annotation(
      Line(points = {{-126, -70}, {-68, -70}}, color = {0, 127, 255}));
    connect(reservoir.ports[3], source.ports[1]) annotation(
      Line(points = {{-126, -70}, {-152, -70}, {-152, -56}, {-182, -56}}, color = {0, 127, 255}));
    connect(pipe.port_b, userValve.port_a) annotation(
      Line(points = {{15, -51}, {15, -30}, {58, -30}}, color = {0, 127, 255}));
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
  end PumpingSystem_reservoirbefore;
  
  model PumpingSystem_co2 "Model of a pumping system for drinking water"
    extends Modelica.Icons.Example;
    replaceable package Medium = ExternalMedia.Examples.CO2CoolProp;
    //replaceable package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater
    //  constrainedby Modelica.Media.Interfaces.PartialMedium;
    Modelica.Fluid.Sources.FixedBoundary source(nPorts = 1, use_T = true, T = Modelica.Units.Conversions.from_degC(20), p (displayUnit = "Pa")= 600000, redeclare package Medium = Medium) annotation(
      Placement(transformation(origin = {-64, 6}, extent = {{-100, -80}, {-80, -60}})));
    Modelica.Fluid.Pipes.StaticPipe pipe(allowFlowReversal = true, length = 100, height_ab = 50, diameter = 0.3, redeclare package Medium = Medium) annotation(
      Placement(transformation(origin = {-30, -51}, extent = {{-9, -10}, {11, 10}}, rotation = 90)));
    Modelica.Fluid.Machines.PrescribedPump pumps(checkValve = true, checkValveHomotopy = Modelica.Fluid.Types.CheckValveHomotopyType.Closed, N_nominal = 1200, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), use_N_in = true, nParallel = 1, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, V(displayUnit = "l") = 0.05, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, redeclare package Medium = Medium, p_b_start = 600000, T_start = system.T_start) annotation(
      Placement(transformation(extent = {{-68, -80}, {-48, -60}})));
    Modelica.Fluid.Valves.ValveLinear userValve(allowFlowReversal = false, dp_nominal = 200000, m_flow_nominal = 400, redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{58, -38}, {74, -22}})));
    Modelica.Fluid.Sources.FixedBoundary sink(p (displayUnit = "Pa")= 600000, T = system.T_ambient, nPorts = 2, redeclare package Medium = Medium) annotation(
      Placement(transformation(extent = {{100, -40}, {80, -20}})));
    inner Modelica.Fluid.System system(energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial) annotation(
      Placement(transformation(extent = {{60, -96}, {80, -76}})));
    pumpSpeedInput pumpSpeedInput1(n_speed = 1200) annotation(
      Placement(transformation(origin = {-56, -44}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    pumpSpeedInput valveopeningconstant(n_speed = 0.5) annotation(
      Placement(transformation(origin = {66, -6}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Fluid.Vessels.OpenTank reservoir(redeclare package Medium = Medium, T_start = Modelica.Units.Conversions.from_degC(20), crossArea = 50, height = 3, level_start = 2.2, nPorts = 3, portsData = {Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.3), Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.3), Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.01)}, use_portsData = true, p_ambient = 6e10) annotation(
      Placement(transformation(origin = {-110, -14}, extent = {{-20, -16}, {0, 4}})));
  equation
    connect(userValve.port_b, sink.ports[1]) annotation(
      Line(points = {{74, -30}, {80, -30}}, color = {0, 127, 255}));
    connect(pipe.port_a, pumps.port_b) annotation(
      Line(points = {{-30, -60}, {-30, -70}, {-48, -70}}, color = {0, 127, 255}));
    connect(pumpSpeedInput1.y, pumps.N_in) annotation(
      Line(points = {{-57, -45}, {-57, -60}, {-58, -60}}, color = {0, 0, 127}));
    connect(valveopeningconstant.y, userValve.opening) annotation(
      Line(points = {{67, -7}, {67, -16}, {72, -16}, {72, -15}, {66, -15}, {66, -24}}, color = {0, 0, 127}));
    connect(reservoir.ports[2], pumps.port_a) annotation(
      Line(points = {{-120, -30}, {-118, -30}, {-118, -36}, {-120, -36}, {-120, -70}, {-68, -70}}, color = {0, 127, 255}));
    connect(reservoir.ports[3], source.ports[1]) annotation(
      Line(points = {{-120, -30}, {-123, -30}, {-123, -64}, {-144, -64}}, color = {0, 127, 255}));
    connect(pipe.port_b, userValve.port_a) annotation(
      Line(points = {{-30, -40}, {58, -40}, {58, -30}}, color = {0, 127, 255}));
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
  end PumpingSystem_co2;
  annotation(
    uses(Modelica(version = "4.0.0")));
end LearningUseCases;
