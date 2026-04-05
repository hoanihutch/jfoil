@testset "panel_linvortex_velocity" begin
    Xj = [0.0 1.0; 0.0 0.0]
    xi = [0.5, 0.1]

    # No vdir, not onmid
    a, b = JFoil.panel_linvortex_velocity(Xj, xi, nothing, false)
    @test a[1] ≈ 0.21858352 atol=1e-6
    @test a[2] ≈ -0.11543824 atol=1e-6
    @test b[1] ≈ 0.21858352 atol=1e-6
    @test b[2] ≈ 0.11543824 atol=1e-6

    # onmid
    xi_mid = [0.5, 0.0]
    a2, b2 = JFoil.panel_linvortex_velocity(Xj, xi_mid, nothing, true)
    @test a2[1] ≈ 0.25 atol=1e-10
    @test a2[2] ≈ -0.15915494 atol=1e-6
    @test b2[1] ≈ 0.25 atol=1e-10
    @test b2[2] ≈ 0.15915494 atol=1e-6

    # with vdir
    vdir = [0.0, 1.0]
    a3, b3 = JFoil.panel_linvortex_velocity(Xj, xi, vdir, false)
    @test a3 ≈ -0.11543823891079545 atol=1e-10
    @test b3 ≈ 0.11543823891079545 atol=1e-10
end

@testset "panel_linvortex_stream" begin
    Xj = [0.0 1.0; 0.0 0.0]
    xi = [0.5, 0.1]
    a, b = JFoil.panel_linvortex_stream(Xj, xi)
    @test a ≈ -0.11131747690107713 atol=1e-10
    @test b ≈ -0.11131747690107713 atol=1e-10
end

@testset "panel_constsource_velocity" begin
    Xj = [0.0 1.0; 0.0 0.0]
    xi = [0.5, 0.1]
    a = JFoil.panel_constsource_velocity(Xj, xi, nothing)
    @test a[1] ≈ 0.0 atol=1e-10
    @test a[2] ≈ 0.43716704 atol=1e-6
end

@testset "panel_constsource_stream" begin
    Xj = [0.0 1.0; 0.0 0.0]
    xi = [0.5, 0.1]
    a = JFoil.panel_constsource_stream(Xj, xi)
    @test a ≈ 0.0 atol=1e-10
end

@testset "panel_linsource_velocity" begin
    Xj = [0.0 1.0; 0.0 0.0]
    xi = [0.5, 0.1]
    a, b = JFoil.panel_linsource_velocity(Xj, xi, nothing)
    @test a[1] ≈ 0.11543824 atol=1e-6
    @test a[2] ≈ 0.21858352 atol=1e-6
    @test b[1] ≈ -0.11543824 atol=1e-6
    @test b[2] ≈ 0.21858352 atol=1e-6
end

@testset "panel_linsource_stream" begin
    Xj = [0.0 1.0; 0.0 0.0]
    xi = [0.5, 0.1]
    a, b = JFoil.panel_linsource_stream(Xj, xi)
    @test a ≈ 0.07612603171916493 atol=1e-10
    @test b ≈ 0.17387396828083512 atol=1e-10
end
