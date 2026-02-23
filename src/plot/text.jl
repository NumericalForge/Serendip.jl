# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

# Cache is required to avoid repeated font-discovery lookups.
const _resolved_font_family_cache = Dict{String, String}()

# Ordered fallback candidates for font-family resolution.
const _font_fallback_candidates = [
    "NewComputerModern",
    "NewComputerModern Regular",
    "Times New Roman",
    "Cambria",
    "Palatino Linotype",
    "Georgia",
    "DejaVu Serif",
    "Liberation Serif",
    "Noto Serif",
    "serif",
    "sans-serif",
]

function _font_family_from_name(name::AbstractString)
    isempty(name) && return nothing
    font = FreeTypeAbstraction.findfont(name)
    font === nothing && return nothing
    return replace(font.family_name, "Math" => "")
end

function resolve_font_family(name::AbstractString)
    key = String(name)
    haskey(_resolved_font_family_cache, key) && return _resolved_font_family_cache[key]

    # First, attempt the requested family.
    family = _font_family_from_name(key)

    # Then try known fallbacks.
    if family === nothing
        for candidate in _font_fallback_candidates
            family = _font_family_from_name(candidate)
            family !== nothing && break
        end
    end

    # Last resort: let Cairo/Pango attempt generic families.
    family = something(family, "serif")
    _resolved_font_family_cache[key] = family

    return family
end

# Backward-compatible alias used across plotting files.
get_font(name::AbstractString) = resolve_font_family(name)

function _current_font_family(cc::CairoContext)
    font_face = ccall((:cairo_get_font_face, Cairo.libcairo), Ptr{Cvoid}, (Ptr{Cvoid},), cc.ptr)
    font_family = ccall((:cairo_toy_font_face_get_family, Cairo.libcairo), Cstring, (Ptr{Cvoid},), font_face)
    family = unsafe_string(font_family)
    return isempty(family) ? resolve_font_family("NewComputerModern") : family
end

function _text_anchor_shift(layout::TypesetLayout, halign::AbstractString, valign::AbstractString)
    rw = halign == "center" ? 0.5 : halign == "right" ? 1.0 : 0.0
    rh = valign == "center" ? 0.5 : valign == "top" ? 1.0 : 0.0

    dx = -rw * layout.width
    dy = -((1 - rh) * layout.max_y + rh * layout.min_y)
    return dx, dy
end

function getsize(cc::CairoContext, str::AbstractString, fontsize::Float64)
    font = _current_font_family(cc)
    layout = layout_typeset(cc, str, fontsize, font=font)
    return layout.width, layout.height
end

function getsize(str::AbstractString, fontsize::Float64)
    surf = CairoImageSurface(4, 4, Cairo.FORMAT_ARGB32)
    cc = CairoContext(surf)
    select_font_face(cc, resolve_font_family("NewComputerModern"), Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
    set_font_size(cc, fontsize)
    return getsize(cc, str, fontsize)
end

function draw_text(cc::CairoContext, x, y, str::AbstractString; halign="center", valign="center", angle=0)
    Cairo.save(cc)

    font = _current_font_family(cc)
    fmatrix = get_font_matrix(cc)
    fontsize = norm([fmatrix.xx, fmatrix.xy])

    layout = layout_typeset(cc, str, fontsize, font=font)
    dx, dy = _text_anchor_shift(layout, halign, valign)

    Cairo.translate(cc, x, y)
    Cairo.rotate(cc, -angle * pi / 180)
    Cairo.translate(cc, dx, dy)

    draw_typeset_layout!(cc, layout, font=font)

    set_font_size(cc, fontsize)
    Cairo.restore(cc)
end
