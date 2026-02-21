@testset "Components" begin
    F = TurboMachineModel.Component
    C = TurboMachineModel.Components
    P = TurboMachineModel.Physics.Fluids
    port_names(c) = Set(p.name for p in F.ports(c))
    has_mapping(c, port_name, label) =
        any(v -> any(m -> m.port == port_name && m.label == label, v.mappings), F.variables(c, :steady))
    mapped_variable_ids(c, sim_type, port_name, label) = Set(
        v.id for v in F.variables(c, sim_type)
        if any(m -> m.port == port_name && m.label == label, v.mappings)
    )
    variable_ids(c, sim_type) = Set(v.id for v in F.variables(c, sim_type))

    cmb = C.Combustor(0.04, 43e6, NamedTuple())
    @test :inlet in port_names(cmb)
    @test :outlet in port_names(cmb)
    @test only(mapped_variable_ids(cmb, :steady, :inlet, :mdot)) ==
          only(mapped_variable_ids(cmb, :steady, :outlet, :mdot))
    @test only(mapped_variable_ids(cmb, :transient, :inlet, :mdot)) ==
          only(mapped_variable_ids(cmb, :transient, :outlet, :mdot))
    @test only(mapped_variable_ids(cmb, :steady, :inlet, :composition)) ==
          only(mapped_variable_ids(cmb, :steady, :outlet, :composition))

    pln = C.Plenum(volume=0.25)
    @test :inlet in port_names(pln)
    @test :outlet in port_names(pln)
    @test has_mapping(pln, :inlet, :pt)
    @test has_mapping(pln, :outlet, :ht)
    @test only(mapped_variable_ids(pln, :steady, :inlet, :pt)) ==
          only(mapped_variable_ids(pln, :steady, :outlet, :pt))
    @test only(mapped_variable_ids(pln, :transient, :inlet, :pt)) ==
          only(mapped_variable_ids(pln, :transient, :outlet, :pt))
    @test only(mapped_variable_ids(pln, :steady, :inlet, :ht)) ==
          only(mapped_variable_ids(pln, :steady, :outlet, :ht))
    @test only(mapped_variable_ids(pln, :transient, :inlet, :ht)) ==
          only(mapped_variable_ids(pln, :transient, :outlet, :ht))
    @test only(mapped_variable_ids(pln, :steady, :inlet, :mdot)) ==
          only(mapped_variable_ids(pln, :steady, :outlet, :mdot))
    @test only(mapped_variable_ids(pln, :transient, :inlet, :mdot)) !=
          only(mapped_variable_ids(pln, :transient, :outlet, :mdot))
    @test only(mapped_variable_ids(pln, :steady, :inlet, :composition)) ==
          only(mapped_variable_ids(pln, :steady, :outlet, :composition))
    @test :retained_mass ∉ variable_ids(pln, :steady)
    @test :retained_mass ∈ variable_ids(pln, :transient)
    @test_throws ErrorException C.Plenum(volume=0.0)

    shaft = C.InertialShaft(J=0.35, damping=0.01, n_ports=3)
    shaft_ports = port_names(shaft)
    @test length(collect(shaft_ports)) == 3
    @test :shaft1 in shaft_ports
    @test :shaft2 in shaft_ports
    @test :shaft3 in shaft_ports
    @test only(mapped_variable_ids(shaft, :steady, :shaft1, :omega)) ==
          only(mapped_variable_ids(shaft, :steady, :shaft2, :omega))
    @test only(mapped_variable_ids(shaft, :steady, :shaft1, :tau)) !=
          only(mapped_variable_ids(shaft, :steady, :shaft2, :tau))

    gb = C.Gearbox(ratio=2.5, efficiency=0.98)
    gb_ports = port_names(gb)
    @test :input in gb_ports
    @test :output in gb_ports
    @test has_mapping(gb, :input, :omega)
    @test has_mapping(gb, :output, :tau)
    @test_throws ErrorException C.Gearbox(ratio=0.0)
    @test_throws ErrorException C.Gearbox(ratio=2.0, efficiency=0.0)

    pm = P.PerformanceMap(
        300.0,
        100_000.0,
        [1.0, 2.0],
        [10.0, 20.0],
        [2.0 3.0; 4.0 5.0],
        [0.8 0.82; 0.9 0.92],
    )
    tm = C.Turbomachine(
        mode=:compressor,
        performance_map=pm,
        eta_guess=0.9,
    )
    @test tm.mode == :compressor
    @test :inlet in port_names(tm)
    @test has_mapping(tm, :inlet, :pt)
    @test only(mapped_variable_ids(tm, :steady, :inlet, :mdot)) ==
          only(mapped_variable_ids(tm, :steady, :outlet, :mdot))
    @test only(mapped_variable_ids(tm, :transient, :inlet, :mdot)) !=
          only(mapped_variable_ids(tm, :transient, :outlet, :mdot))
end
