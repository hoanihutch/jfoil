@testset "norm2" begin
    @test JFoil.norm2([3.0, 4.0]) ≈ 5.0
    @test JFoil.norm2([1.0, 0.0]) ≈ 1.0
    @test JFoil.norm2([0.5, -0.5]) ≈ 0.7071067811865476
    @test JFoil.norm2([0.0, 0.0]) ≈ 0.0
end

@testset "dist" begin
    @test JFoil.dist(3.0, 4.0) ≈ 5.0
    @test JFoil.dist(0.0, 0.0) ≈ 0.0
    @test JFoil.dist(1.0, 1.0) ≈ 1.4142135623730951
end

@testset "vprint" begin
    param = JFoil.Param()
    # Default verb=1, so verb=1 should print, verb=2 should not
    mktemp() do path, io
        redirect_stdout(io) do
            JFoil.vprint(param, 1, "test message")
        end
        flush(io)
        @test occursin("test message", read(path, String))
    end

    mktemp() do path, io
        redirect_stdout(io) do
            JFoil.vprint(param, 2, "should not print")
        end
        flush(io)
        @test read(path, String) == ""
    end
end

@testset "Struct defaults" begin
    g = JFoil.Geom()
    @test g.chord == 1.0
    @test g.wakelen == 1.0
    @test g.npoint == 1
    @test g.name == "noname"
    @test g.xref == [0.25, 0.0]

    p = JFoil.Panel()
    @test p.N == 0

    o = JFoil.Oper()
    @test o.Vinf == 1.0
    @test o.alpha == 0.0
    @test o.Re == 1e5
    @test o.Ma == 0.0
    @test o.viscous == false

    i = JFoil.Isol()
    @test i.sstag == 0.0
    @test i.sstag_g == [0.0, 0.0]

    v = JFoil.Vsol()
    @test v.xt == 0.0
    @test size(v.Xt) == (2, 2)

    gl = JFoil.Glob()
    @test gl.Nsys == 0
    @test gl.conv == true

    po = JFoil.Post()
    @test po.cl == 0.0
    @test po.cd == 0.0

    pa = JFoil.Param()
    @test pa.verb == 1
    @test pa.rtol == 1e-10
    @test pa.niglob == 50
    @test pa.ncrit == 9.0
    @test pa.GA == 6.7
    @test pa.GB == 0.75
    @test pa.GC == 18.0
    @test pa.gam == 1.4
    @test pa.Tsrat == 0.35

    m = JFoil.Mfoil()
    @test m.version == "2022-02-22"
    @test m.geom isa JFoil.Geom
    @test m.param.ncrit == 9.0
end
