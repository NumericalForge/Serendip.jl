using Serendip
using Test


struct XmlTestRecord
    label::Any
    values::Vector{Int}
end


@announced_testset "Lightweight XML" begin
    @announced_testset "Parsing" begin
        xml = "﻿<?xml version='1.0' encoding=\"UTF-8\"?>\n" *
              "<!--generated > safely-->\n" *
              "<model-root xmlns:serendip=\"urn:serendip\" enabled = 'true' ν='0.2'>\n" *
              "  <serendip:item id=\"1\">A &amp; B &lt; C &#x1f642;</serendip:item>\n" *
              "  <σ-node/>\n" *
              "</model-root>"

        doc = XmlDocument(xml)
        @test doc.attributes == Dict("version" => "1.0", "encoding" => "UTF-8")
        @test length(doc.children) == 2
        @test doc.root.name == "model-root"
        @test doc.root.attributes["xmlns:serendip"] == "urn:serendip"
        @test doc.root.attributes["enabled"] == "true"
        @test doc.root.attributes["ν"] == "0.2"
        @test doc.root("serendip:item").content == "A & B < C 🙂"
        @test doc.root("σ-node") isa XmlElement
        @test doc.root["id" => "1"] == [doc.root("serendip:item")]

        no_prolog = XmlDocument("  <root><child/></root>")
        @test isempty(no_prolog.attributes)
        @test no_prolog.root.name == "root"

        mktempdir() do dir
            filename = joinpath(dir, "no-prolog.xml")
            save(no_prolog, filename)
            @test !startswith(read(filename, String), "<?xml")
            @test XmlDocument(filename).root("child") isa XmlElement
        end
    end

    @announced_testset "Writing and round-trip" begin
        mktempdir() do dir
            root = XmlElement(
                "model-root",
                attributes=Dict("quoted" => "A&B<\"'", "ν" => "0.2"),
                children=XmlElement[
                    XmlElement("message", content="ação & cálculo < 1"),
                    XmlElement("empty-node"),
                ],
            )
            doc = XmlDocument(("version" => "1.0", "encoding" => "UTF-8"))
            Serendip.addchild!(doc, Serendip.XmlComment("generated > safely"))
            Serendip.addchild!(doc, root)

            filename = joinpath(dir, "document.xml")
            @test save(doc, filename) === nothing
            saved = read(filename, String)
            @test occursin("A&amp;B&lt;&quot;&apos;", saved)
            @test occursin("ação &amp; cálculo &lt; 1", saved)

            loaded = XmlDocument(filename)
            @test length(loaded.children) == 2
            @test loaded.root.attributes["quoted"] == "A&B<\"'"
            @test loaded.root("message").content == "ação & cálculo < 1"

            payload = UInt8[0x00, 0xff, UInt8('<'), UInt8('&'), 0x0a]
            raw = XmlElement("AppendedData", attributes=Dict("encoding" => "raw"), content=payload)
            raw_doc = XmlDocument(("version" => "1.0",), raw)
            raw_filename = joinpath(dir, "raw.xml")
            save(raw_doc, raw_filename)
            @test XmlDocument(raw_filename).root.content == payload
        end
    end

    @announced_testset "Invalid or unsupported input" begin
        invalid = [
            "<first/><second/>",
            "<root><child></root>",
            "<root a='1' a='2'/>",
            "<root>&unknown;</root>",
            "<root>before<child/>after</root>",
            "<root><![CDATA[text]]></root>",
            "<!DOCTYPE root><root/>",
            "<?other value='1'?><root/>",
            "<?xml encoding='UTF-8'?><root/>",
            "<?xml version='1.0' encoding='ISO-8859-1'?><root/>",
            "<!-- unclosed<root/>",
            "<root>invalid\x01control</root>",
        ]
        for xml in invalid
            @test_throws ErrorException XmlDocument(xml)
        end

        mktempdir() do dir
            invalid_name = XmlDocument(("version" => "1.0",), XmlElement("invalid name"))
            @test_throws ErrorException save(invalid_name, joinpath(dir, "invalid-name.xml"))

            invalid_comment = XmlDocument(("version" => "1.0",))
            Serendip.addchild!(invalid_comment, Serendip.XmlComment("not--valid"))
            Serendip.addchild!(invalid_comment, XmlElement("root"))
            @test_throws ErrorException save(invalid_comment, joinpath(dir, "invalid-comment.xml"))

            mixed = XmlElement("root", children=XmlElement[XmlElement("child")], content="text")
            @test_throws ErrorException save(
                XmlDocument(("version" => "1.0",), mixed),
                joinpath(dir, "mixed.xml"),
            )
        end
    end

    @announced_testset "to_xml_node" begin
        attributes = Dict("source" => "test")
        array_node = to_xml_node(view([1, 2, 3], :), "values", attributes)
        @test array_node.attributes["format"] == "ascii"
        @test array_node.attributes["components"] == "1"
        @test array_node.content == "1\n2\n3"
        @test attributes == Dict("source" => "test")

        heterogeneous = to_xml_node(Any[1, "two"])
        @test length(heterogeneous.children) == 2
        @test heterogeneous.children[1].content == "1"
        @test heterogeneous.children[2].content == "two"

        record = to_xml_node(XmlTestRecord("sample", [2, 4]))
        @test record.attributes["label"] == "sample"
        @test record("values").content == "2\n4"
    end
end
