using Serendip

geo = GeoModel()
add_block(geo, [0.0, 0.0, 0.0], 1, 1, 1, nx=2, ny=2, nz=2, shape=HEX8)

mesh = Mesh(geo)

save(mesh, "mesh0.vtu")

select(mesh, :element, x>0.25, tag="right")

# Serendip.add_contact_elements(mesh)
# save(mesh, "mesh1.vtu")

Serendip.add_cohesive_elements(mesh)
save(mesh, "mesh1.vtu")