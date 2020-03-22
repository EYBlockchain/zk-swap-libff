#include <libff/algebra/curves/bw12_446/bw12_446_g2.hpp>

namespace libff {

#ifdef PROFILE_OP_COUNTS
long long bw12_446_G2::add_cnt = 0;
long long bw12_446_G2::dbl_cnt = 0;
#endif

std::vector<size_t> bw12_446_G2::wnaf_window_table;
std::vector<size_t> bw12_446_G2::fixed_base_exp_window_table;
bw12_446_G2 bw12_446_G2::G2_zero;
bw12_446_G2 bw12_446_G2::G2_one;

bw12_446_G2::bw12_446_G2()
{
    this->X = G2_zero.X;
    this->Y = G2_zero.Y;
    this->Z = G2_zero.Z;
}

bw12_446_Fq2 bw12_446_G2::mul_by_b(const bw12_446_Fq2 &elt)
{
    return bw12_446_Fq2(bw12_446_twist_mul_by_b_c0 * elt.c0, bw12_446_twist_mul_by_b_c1 * elt.c1);
}

void bw12_446_G2::print() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        bw12_446_G2 copy(*this);
        copy.to_affine_coordinates();
        gmp_printf("(%Nd*z + %Nd , %Nd*z + %Nd)\n",
                   copy.X.c1.as_bigint().data, bw12_446_Fq::num_limbs,
                   copy.X.c0.as_bigint().data, bw12_446_Fq::num_limbs,
                   copy.Y.c1.as_bigint().data, bw12_446_Fq::num_limbs,
                   copy.Y.c0.as_bigint().data, bw12_446_Fq::num_limbs);
    }
}

void bw12_446_G2::print_coordinates() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        gmp_printf("(%Nd*z + %Nd : %Nd*z + %Nd : %Nd*z + %Nd)\n",
                   this->X.c1.as_bigint().data, bw12_446_Fq::num_limbs,
                   this->X.c0.as_bigint().data, bw12_446_Fq::num_limbs,
                   this->Y.c1.as_bigint().data, bw12_446_Fq::num_limbs,
                   this->Y.c0.as_bigint().data, bw12_446_Fq::num_limbs,
                   this->Z.c1.as_bigint().data, bw12_446_Fq::num_limbs,
                   this->Z.c0.as_bigint().data, bw12_446_Fq::num_limbs);
    }
}

void bw12_446_G2::to_affine_coordinates()
{
    if (this->is_zero())
    {
        this->X = bw12_446_Fq2::zero();
        this->Y = bw12_446_Fq2::one();
        this->Z = bw12_446_Fq2::zero();
    }
    else
    {
        bw12_446_Fq2 Z_inv = Z.inverse();
        bw12_446_Fq2 Z2_inv = Z_inv.squared();
        bw12_446_Fq2 Z3_inv = Z2_inv * Z_inv;
        this->X = this->X * Z2_inv;
        this->Y = this->Y * Z3_inv;
        this->Z = bw12_446_Fq2::one();
    }
}

void bw12_446_G2::to_special()
{
    this->to_affine_coordinates();
}

bool bw12_446_G2::is_special() const
{
    return (this->is_zero() || this->Z == bw12_446_Fq2::one());
}

bool bw12_446_G2::is_zero() const
{
    return (this->Z.is_zero());
}

bool bw12_446_G2::operator==(const bw12_446_G2 &other) const
{
    if (this->is_zero())
    {
        return other.is_zero();
    }

    if (other.is_zero())
    {
        return false;
    }

    /* now neither is O */

    // using Jacobian coordinates so:
    // (X1:Y1:Z1) = (X2:Y2:Z2)
    // iff
    // X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
    // iff
    // X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3

    bw12_446_Fq2 Z1_squared = (this->Z).squared();
    bw12_446_Fq2 Z2_squared = (other.Z).squared();

    if ((this->X * Z2_squared) != (other.X * Z1_squared))
    {
        return false;
    }

    bw12_446_Fq2 Z1_cubed = (this->Z) * Z1_squared;
    bw12_446_Fq2 Z2_cubed = (other.Z) * Z2_squared;

    if ((this->Y * Z2_cubed) != (other.Y * Z1_cubed))
    {
        return false;
    }

    return true;
}

bool bw12_446_G2::operator!=(const bw12_446_G2& other) const
{
    return !(operator==(other));
}

bw12_446_G2 bw12_446_G2::operator+(const bw12_446_G2 &other) const
{
    // handle special cases having to do with O
    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return *this;
    }

    // no need to handle points of order 2,4
    // (they cannot exist in a prime-order subgroup)

    // check for doubling case

    // using Jacobian coordinates so:
    // (X1:Y1:Z1) = (X2:Y2:Z2)
    // iff
    // X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
    // iff
    // X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3

    bw12_446_Fq2 Z1Z1 = (this->Z).squared();
    bw12_446_Fq2 Z2Z2 = (other.Z).squared();

    bw12_446_Fq2 U1 = this->X * Z2Z2;
    bw12_446_Fq2 U2 = other.X * Z1Z1;

    bw12_446_Fq2 Z1_cubed = (this->Z) * Z1Z1;
    bw12_446_Fq2 Z2_cubed = (other.Z) * Z2Z2;

    bw12_446_Fq2 S1 = (this->Y) * Z2_cubed;      // S1 = Y1 * Z2 * Z2Z2
    bw12_446_Fq2 S2 = (other.Y) * Z1_cubed;      // S2 = Y2 * Z1 * Z1Z1

    if (U1 == U2 && S1 == S2)
    {
        // dbl case; nothing of above can be reused
        return this->dbl();
    }

    // rest of add case
    bw12_446_Fq2 H = U2 - U1;                            // H = U2-U1
    bw12_446_Fq2 S2_minus_S1 = S2-S1;
    bw12_446_Fq2 I = (H+H).squared();                    // I = (2 * H)^2
    bw12_446_Fq2 J = H * I;                              // J = H * I
    bw12_446_Fq2 r = S2_minus_S1 + S2_minus_S1;          // r = 2 * (S2-S1)
    bw12_446_Fq2 V = U1 * I;                             // V = U1 * I
    bw12_446_Fq2 X3 = r.squared() - J - (V+V);           // X3 = r^2 - J - 2 * V
    bw12_446_Fq2 S1_J = S1 * J;
    bw12_446_Fq2 Y3 = r * (V-X3) - (S1_J+S1_J);          // Y3 = r * (V-X3)-2 S1 J
    bw12_446_Fq2 Z3 = ((this->Z+other.Z).squared()-Z1Z1-Z2Z2) * H; // Z3 = ((Z1+Z2)^2-Z1Z1-Z2Z2) * H

    return bw12_446_G2(X3, Y3, Z3);
}

bw12_446_G2 bw12_446_G2::operator-() const
{
    return bw12_446_G2(this->X, -(this->Y), this->Z);
}


bw12_446_G2 bw12_446_G2::operator-(const bw12_446_G2 &other) const
{
    return (*this) + (-other);
}

bw12_446_G2 bw12_446_G2::add(const bw12_446_G2 &other) const
{
    // handle special cases having to do with O
    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return *this;
    }

    // no need to handle points of order 2,4
    // (they cannot exist in a prime-order subgroup)

    // handle double case
    if (this->operator==(other))
    {
        return this->dbl();
    }

#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif
    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#addition-add-1998-cmo-2

    bw12_446_Fq2 Z1Z1 = (this->Z).squared();             // Z1Z1 = Z1^2
    bw12_446_Fq2 Z2Z2 = (other.Z).squared();             // Z2Z2 = Z2^2
    bw12_446_Fq2 U1 = (this->X) * Z2Z2;                  // U1 = X1 * Z2Z2
    bw12_446_Fq2 U2 = (other.X) * Z1Z1;                  // U2 = X2 * Z1Z1
    bw12_446_Fq2 S1 = (this->Y) * (other.Z) * Z2Z2;      // S1 = Y1 * Z2 * Z2Z2
    bw12_446_Fq2 S2 = (other.Y) * (this->Z) * Z1Z1;      // S2 = Y2 * Z1 * Z1Z1
    bw12_446_Fq2 H = U2 - U1;                            // H = U2-U1
    bw12_446_Fq2 S2_minus_S1 = S2-S1;
    bw12_446_Fq2 I = (H+H).squared();                    // I = (2 * H)^2
    bw12_446_Fq2 J = H * I;                              // J = H * I
    bw12_446_Fq2 r = S2_minus_S1 + S2_minus_S1;          // r = 2 * (S2-S1)
    bw12_446_Fq2 V = U1 * I;                             // V = U1 * I
    bw12_446_Fq2 X3 = r.squared() - J - (V+V);           // X3 = r^2 - J - 2 * V
    bw12_446_Fq2 S1_J = S1 * J;
    bw12_446_Fq2 Y3 = r * (V-X3) - (S1_J+S1_J);          // Y3 = r * (V-X3)-2 S1 J
    bw12_446_Fq2 Z3 = ((this->Z+other.Z).squared()-Z1Z1-Z2Z2) * H; // Z3 = ((Z1+Z2)^2-Z1Z1-Z2Z2) * H

    return bw12_446_G2(X3, Y3, Z3);
}

bw12_446_G2 bw12_446_G2::mixed_add(const bw12_446_G2 &other) const
{
#ifdef DEBUG
    assert(other.is_special());
#endif

    // handle special cases having to do with O
    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return *this;
    }

    // no need to handle points of order 2,4
    // (they cannot exist in a prime-order subgroup)

    // check for doubling case

    // using Jacobian coordinates so:
    // (X1:Y1:Z1) = (X2:Y2:Z2)
    // iff
    // X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
    // iff
    // X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3

    // we know that Z2 = 1

    const bw12_446_Fq2 Z1Z1 = (this->Z).squared();

    const bw12_446_Fq2 &U1 = this->X;
    const bw12_446_Fq2 U2 = other.X * Z1Z1;

    const bw12_446_Fq2 Z1_cubed = (this->Z) * Z1Z1;

    const bw12_446_Fq2 &S1 = (this->Y);                // S1 = Y1 * Z2 * Z2Z2
    const bw12_446_Fq2 S2 = (other.Y) * Z1_cubed;      // S2 = Y2 * Z1 * Z1Z1

    if (U1 == U2 && S1 == S2)
    {
        // dbl case; nothing of above can be reused
        return this->dbl();
    }

#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif

    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-madd-2007-bl
    bw12_446_Fq2 H = U2-(this->X);                         // H = U2-X1
    bw12_446_Fq2 HH = H.squared() ;                        // HH = H&2
    bw12_446_Fq2 I = HH+HH;                                // I = 4*HH
    I = I + I;
    bw12_446_Fq2 J = H*I;                                  // J = H*I
    bw12_446_Fq2 r = S2-(this->Y);                         // r = 2*(S2-Y1)
    r = r + r;
    bw12_446_Fq2 V = (this->X) * I ;                       // V = X1*I
    bw12_446_Fq2 X3 = r.squared()-J-V-V;                   // X3 = r^2-J-2*V
    bw12_446_Fq2 Y3 = (this->Y)*J;                         // Y3 = r*(V-X3)-2*Y1*J
    Y3 = r*(V-X3) - Y3 - Y3;
    bw12_446_Fq2 Z3 = ((this->Z)+H).squared() - Z1Z1 - HH; // Z3 = (Z1+H)^2-Z1Z1-HH

    return bw12_446_G2(X3, Y3, Z3);
}

bw12_446_G2 bw12_446_G2::dbl() const
{
#ifdef PROFILE_OP_COUNTS
    this->dbl_cnt++;
#endif
    // handle point at infinity
    if (this->is_zero())
    {
        return (*this);
    }

    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#doubling-dbl-2007-bl

    bw12_446_Fq2 A = (this->X).squared();         // A = X1^2
    bw12_446_Fq2 B = (this->Y).squared();        // B = Y1^2
    bw12_446_Fq2 C = B.squared();                // C = B^2
    bw12_446_Fq2 D = (this->X + B).squared() - A - C;
    D = D+D;                        // D = 2 * ((X1 + B)^2 - A - C)
    bw12_446_Fq2 E = A + A + A;                  // E = 3 * A
    bw12_446_Fq2 F = E.squared();                // F = E^2
    bw12_446_Fq2 X3 = F - (D+D);                 // X3 = F - 2 D
    bw12_446_Fq2 eightC = C+C;
    eightC = eightC + eightC;
    eightC = eightC + eightC;
    bw12_446_Fq2 Y3 = E * (D - X3) - eightC;     // Y3 = E * (D - X3) - 8 * C
    bw12_446_Fq2 Y1Z1 = (this->Y)*(this->Z);
    bw12_446_Fq2 Z3 = Y1Z1 + Y1Z1;               // Z3 = 2 * Y1 * Z1

    return bw12_446_G2(X3, Y3, Z3);
}

// bw12_446_G2 bw12_446_G2::mul_by_q() const
// {
//     return bw12_446_G2::base_field_char() * (*this);
// }


bw12_446_G2 bw12_446_G2::mul_by_q() const
{
    return bw12_446_G2(bw12_446_twist_mul_by_q_X * (this->X).Frobenius_map(1),
                      bw12_446_twist_mul_by_q_Y * (this->Y).Frobenius_map(1),
                      (this->Z).Frobenius_map(1));
}

bool bw12_446_G2::is_well_formed() const
{
    if (this->is_zero())
    {
        return true;
    }
    else
    {
        /*
          y^2 = x^3 + b

          We are using Jacobian coordinates, so equation we need to check is actually

          (y/z^3)^2 = (x/z^2)^3 + b
          y^2 / z^6 = x^3 / z^6 + b
          y^2 = x^3 + b z^6
        */
        bw12_446_Fq2 X2 = this->X.squared();
        bw12_446_Fq2 Y2 = this->Y.squared();
        bw12_446_Fq2 Z2 = this->Z.squared();

        bw12_446_Fq2 X3 = this->X * X2;
        bw12_446_Fq2 Z3 = this->Z * Z2;
        bw12_446_Fq2 Z6 = Z3.squared();

        return (Y2 == X3 + bw12_446_twist_coeff_b * Z6);
    }
}

bw12_446_G2 bw12_446_G2::zero()
{
    return G2_zero;
}

bw12_446_G2 bw12_446_G2::one()
{
    return G2_one;
}

bw12_446_G2 bw12_446_G2::random_element()
{
    return (bw12_446_Fr::random_element().as_bigint()) * G2_one;
}

std::ostream& operator<<(std::ostream &out, const bw12_446_G2 &g)
{
    bw12_446_G2 copy(g);
    copy.to_affine_coordinates();
    out << (copy.is_zero() ? 1 : 0) << OUTPUT_SEPARATOR;
#ifdef NO_PT_COMPRESSION
    out << copy.X << OUTPUT_SEPARATOR << copy.Y;
#else
    /* storing LSB of Y */
    out << copy.X << OUTPUT_SEPARATOR << (copy.Y.c0.as_bigint().data[0] & 1);
#endif

    return out;
}

std::istream& operator>>(std::istream &in, bw12_446_G2 &g)
{
    char is_zero;
    bw12_446_Fq2 tX, tY;

#ifdef NO_PT_COMPRESSION
    in >> is_zero >> tX >> tY;
    is_zero -= '0';
#else
    in.read((char*)&is_zero, 1); // this reads is_zero;
    is_zero -= '0';
    consume_OUTPUT_SEPARATOR(in);

    unsigned char Y_lsb;
    in >> tX;
    consume_OUTPUT_SEPARATOR(in);
    in.read((char*)&Y_lsb, 1);
    Y_lsb -= '0';

    // y = +/- sqrt(x^3 + b)
    if (!is_zero)
    {
        bw12_446_Fq2 tX2 = tX.squared();
        bw12_446_Fq2 tY2 = tX2 * tX + bw12_446_twist_coeff_b;
        tY = tY2.sqrt();

        if ((tY.c0.as_bigint().data[0] & 1) != Y_lsb)
        {
            tY = -tY;
        }
    }
#endif
    // using projective coordinates
    if (!is_zero)
    {
        g.X = tX;
        g.Y = tY;
        g.Z = bw12_446_Fq2::one();
    }
    else
    {
        g = bw12_446_G2::zero();
    }

    return in;
}

void bw12_446_G2::batch_to_special_all_non_zeros(std::vector<bw12_446_G2> &vec)
{
    std::vector<bw12_446_Fq2> Z_vec;
    Z_vec.reserve(vec.size());

    for (auto &el: vec)
    {
        Z_vec.emplace_back(el.Z);
    }
    batch_invert<bw12_446_Fq2>(Z_vec);

    const bw12_446_Fq2 one = bw12_446_Fq2::one();

    for (size_t i = 0; i < vec.size(); ++i)
    {
        bw12_446_Fq2 Z2 = Z_vec[i].squared();
        bw12_446_Fq2 Z3 = Z_vec[i] * Z2;

        vec[i].X = vec[i].X * Z2;
        vec[i].Y = vec[i].Y * Z3;
        vec[i].Z = one;
    }
}

} // libff
