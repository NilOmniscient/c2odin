# c2odin
Port of Randy Gaul's cute_headers/cute_c2 collision header. 
Original source here: https://github.com/RandyGaul/cute_headers/blob/master/cute_c2.h
Original commit used: ad9cce337f5ca93196b4f16f6b358cbbf03124ef

Honestly, I just ported this to better teach myself how to work in Odin. 
Code is under the original source licences as it's merely a port.

If anyone actually wants to use this, the important bits 
are the shapes, and the collision functions. 

# Shapes
* Circle
* Capsule
* AABB
* Polygon
* Ray

# Collision Functions
* Collide       -> Constructs a manifold to describe how the shapes hit
* Collided      -> bool yes/no collision check. 
* CastRay       -> Casts a ray at a shape and returns collision details. 
* Create<SHAPE> -> Creates the various shapes
* GJK           -> Runs the GJK algorithm to calculate closest point between two shapes.
* TOI           -> Computes the time of impact between two shapes. 

For other information, please check out the original. 