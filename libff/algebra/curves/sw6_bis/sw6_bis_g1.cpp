#include <libff/algebra/curves/sw6_bis/sw6_bis_g1.hpp>

namespace libff {

#ifdef PROFILE_OP_COUNTS
long long sw6_bis_G1::add_cnt = 0;
long long sw6_bis_G1::dbl_cnt = 0;
#endif

std::vector<size_t> sw6_bis_G1::wnaf_window_table;
std::vector<size_t> sw6_bis_G1::fixed_base_exp_window_table;
sw6_bis_G1 sw6_bis_G1::G1_zero;
sw6_bis_G1 sw6_bis_G1::G1_one;
sw6_bis_Fq sw6_bis_G1::coeff_a;
sw6_bis_Fq sw6_bis_G1::coeff_b;

sw6_bis_G1::sw6_bis_G1()
{
    this->X_ = G1_zero.X_;
    this->Y_ = G1_zero.Y_;
    this->Z_ = G1_zero.Z_;
}

void sw6_bis_G1::print() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        sw6_bis_G1 copy(*this);
        copy.to_affine_coordinates();
        gmp_printf("(%Nd , %Nd)\n",
                   copy.X_.as_bigint().data, sw6_bis_Fq::num_limbs,
                   copy.Y_.as_bigint().data, sw6_bis_Fq::num_limbs);
    }
}

void sw6_bis_G1::print_coordinates() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        gmp_printf("(%Nd : %Nd : %Nd)\n",
                   this->X_.as_bigint().data, sw6_bis_Fq::num_limbs,
                   this->Y_.as_bigint().data, sw6_bis_Fq::num_limbs,
                   this->Z_.as_bigint().data, sw6_bis_Fq::num_limbs);
    }
}

void sw6_bis_G1::to_affine_coordinates()
{
    if (this->is_zero())
    {
        this->X_ = sw6_bis_Fq::zero();
        this->Y_ = sw6_bis_Fq::one();
        this->Z_ = sw6_bis_Fq::zero();
    }
    else
    {
        const sw6_bis_Fq Z_inv = Z_.inverse();
        this->X_ = this->X_ * Z_inv;
        this->Y_ = this->Y_ * Z_inv;
        this->Z_ = sw6_bis_Fq::one();
    }
}

void sw6_bis_G1::to_special()
{
    this->to_affine_coordinates();
}

bool sw6_bis_G1::is_special() const
{
    return (this->is_zero() || this->Z_ == sw6_bis_Fq::one());
}

bool sw6_bis_G1::is_zero() const
{
    return (this->X_.is_zero() && this->Z_.is_zero());
}

bool sw6_bis_G1::operator==(const sw6_bis_G1 &other) const
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

    // X1/Z1 = X2/Z2 <=> X1*Z2 = X2*Z1
    if ((this->X_ * other.Z_) != (other.X_ * this->Z_))
    {
        return false;
    }

    // Y1/Z1 = Y2/Z2 <=> Y1*Z2 = Y2*Z1
    if ((this->Y_ * other.Z_) != (other.Y_ * this->Z_))
    {
        return false;
    }

    return true;
}

bool sw6_bis_G1::operator!=(const sw6_bis_G1& other) const
{
    return !(operator==(other));
}

sw6_bis_G1 sw6_bis_G1::operator+(const sw6_bis_G1 &other) const
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

    // handle double case, and then all the rest
    /*
      The code below is equivalent to (but faster than) the snippet below:

      if (this->operator==(other))
      {
      return this->dbl();
      }
      else
      {
      return this->add(other);
      }
    */

    const sw6_bis_Fq X1Z2 = (this->X_) * (other.Z_);        // X1Z2 = X1*Z2
    const sw6_bis_Fq X2Z1 = (this->Z_) * (other.X_);        // X2Z1 = X2*Z1

    // (used both in add and double checks)

    const sw6_bis_Fq Y1Z2 = (this->Y_) * (other.Z_);        // Y1Z2 = Y1*Z2
    const sw6_bis_Fq Y2Z1 = (this->Z_) * (other.Y_);        // Y2Z1 = Y2*Z1

    if (X1Z2 == X2Z1 && Y1Z2 == Y2Z1)
    {
        // perform dbl case
        const sw6_bis_Fq XX   = (this->X_).squared();                   // XX  = X1^2
        const sw6_bis_Fq ZZ   = (this->Z_).squared();                   // ZZ  = Z1^2
        const sw6_bis_Fq w    = sw6_bis_G1::coeff_a * ZZ + (XX + XX + XX); // w   = a*ZZ + 3*XX
        const sw6_bis_Fq Y1Z1 = (this->Y_) * (this->Z_);
        const sw6_bis_Fq s    = Y1Z1 + Y1Z1;                            // s   = 2*Y1*Z1
        const sw6_bis_Fq ss   = s.squared();                            // ss  = s^2
        const sw6_bis_Fq sss  = s * ss;                                 // sss = s*ss
        const sw6_bis_Fq R    = (this->Y_) * s;                         // R   = Y1*s
        const sw6_bis_Fq RR   = R.squared();                            // RR  = R^2
        const sw6_bis_Fq B    = ((this->X_)+R).squared()-XX-RR;         // B   = (X1+R)^2 - XX - RR
        const sw6_bis_Fq h    = w.squared() - (B+B);                    // h   = w^2 - 2*B
        const sw6_bis_Fq X3   = h * s;                                  // X3  = h*s
        const sw6_bis_Fq Y3   = w * (B-h)-(RR+RR);                      // Y3  = w*(B-h) - 2*RR
        const sw6_bis_Fq Z3   = sss;                                    // Z3  = sss

        return sw6_bis_G1(X3, Y3, Z3);
    }

    // if we have arrived here we are in the add case
    const sw6_bis_Fq Z1Z2 = (this->Z_) * (other.Z_);      // Z1Z2 = Z1*Z2
    const sw6_bis_Fq u    = Y2Z1 - Y1Z2;                  // u    = Y2*Z1-Y1Z2
    const sw6_bis_Fq uu   = u.squared();                  // uu   = u^2
    const sw6_bis_Fq v    = X2Z1 - X1Z2;                  // v    = X2*Z1-X1Z2
    const sw6_bis_Fq vv   = v.squared();                  // vv   = v^2
    const sw6_bis_Fq vvv  = v * vv;                       // vvv  = v*vv
    const sw6_bis_Fq R    = vv * X1Z2;                    // R    = vv*X1Z2
    const sw6_bis_Fq A    = uu * Z1Z2 - (vvv + R + R);    // A    = uu*Z1Z2 - vvv - 2*R
    const sw6_bis_Fq X3   = v * A;                        // X3   = v*A
    const sw6_bis_Fq Y3   = u * (R-A) - vvv * Y1Z2;       // Y3   = u*(R-A) - vvv*Y1Z2
    const sw6_bis_Fq Z3   = vvv * Z1Z2;                   // Z3   = vvv*Z1Z2

    return sw6_bis_G1(X3, Y3, Z3);
}

sw6_bis_G1 sw6_bis_G1::operator-() const
{
    return sw6_bis_G1(this->X_, -(this->Y_), this->Z_);
}


sw6_bis_G1 sw6_bis_G1::operator-(const sw6_bis_G1 &other) const
{
    return (*this) + (-other);
}

sw6_bis_G1 sw6_bis_G1::add(const sw6_bis_G1 &other) const
{
    // handle special cases having to do with O
    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return (*this);
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

    const sw6_bis_Fq Y1Z2 = (this->Y_) * (other.Z_);        // Y1Z2 = Y1*Z2
    const sw6_bis_Fq X1Z2 = (this->X_) * (other.Z_);        // X1Z2 = X1*Z2
    const sw6_bis_Fq Z1Z2 = (this->Z_) * (other.Z_);        // Z1Z2 = Z1*Z2
    const sw6_bis_Fq u    = (other.Y_) * (this->Z_) - Y1Z2; // u    = Y2*Z1-Y1Z2
    const sw6_bis_Fq uu   = u.squared();                    // uu   = u^2
    const sw6_bis_Fq v    = (other.X_) * (this->Z_) - X1Z2; // v    = X2*Z1-X1Z2
    const sw6_bis_Fq vv   = v.squared();                    // vv   = v^2
    const sw6_bis_Fq vvv  = v * vv;                         // vvv  = v*vv
    const sw6_bis_Fq R    = vv * X1Z2;                      // R    = vv*X1Z2
    const sw6_bis_Fq A    = uu * Z1Z2 - (vvv + R + R);      // A    = uu*Z1Z2 - vvv - 2*R
    const sw6_bis_Fq X3   = v * A;                          // X3   = v*A
    const sw6_bis_Fq Y3   = u * (R-A) - vvv * Y1Z2;         // Y3   = u*(R-A) - vvv*Y1Z2
    const sw6_bis_Fq Z3   = vvv * Z1Z2;                     // Z3   = vvv*Z1Z2

    return sw6_bis_G1(X3, Y3, Z3);
}

sw6_bis_G1 sw6_bis_G1::mixed_add(const sw6_bis_G1 &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif
    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#addition-add-1998-cmo-2
    //assert(other.Z == sw6_bis_Fq::one());

    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return (*this);
    }

#ifdef DEBUG
    assert(other.is_special());
#endif

    const sw6_bis_Fq &X1Z2 = (this->X_);                    // X1Z2 = X1*Z2 (but other is special and not zero)
    const sw6_bis_Fq X2Z1 = (this->Z_) * (other.X_);        // X2Z1 = X2*Z1

    // (used both in add and double checks)

    const sw6_bis_Fq &Y1Z2 = (this->Y_);                    // Y1Z2 = Y1*Z2 (but other is special and not zero)
    const sw6_bis_Fq Y2Z1 = (this->Z_) * (other.Y_);        // Y2Z1 = Y2*Z1

    if (X1Z2 == X2Z1 && Y1Z2 == Y2Z1)
    {
        return this->dbl();
    }

    sw6_bis_Fq u = Y2Z1 - this->Y_;             // u = Y2*Z1-Y1
    sw6_bis_Fq uu = u.squared();                // uu = u2
    sw6_bis_Fq v = X2Z1 - this->X_;             // v = X2*Z1-X1
    sw6_bis_Fq vv = v.squared();                // vv = v2
    sw6_bis_Fq vvv = v*vv;                      // vvv = v*vv
    sw6_bis_Fq R = vv * this->X_;               // R = vv*X1
    sw6_bis_Fq A = uu * this->Z_ - vvv - R - R; // A = uu*Z1-vvv-2*R
    sw6_bis_Fq X3 = v * A;                      // X3 = v*A
    sw6_bis_Fq Y3 = u*(R-A) - vvv * this->Y_;   // Y3 = u*(R-A)-vvv*Y1
    sw6_bis_Fq Z3 = vvv * this->Z_;             // Z3 = vvv*Z1

    return sw6_bis_G1(X3, Y3, Z3);
}

sw6_bis_G1 sw6_bis_G1::dbl() const
{
#ifdef PROFILE_OP_COUNTS
    this->dbl_cnt++;
#endif
    if (this->is_zero())
    {
        return (*this);
    }
    else
    {
        // NOTE: does not handle O and pts of order 2,4
        // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#doubling-dbl-2007-bl

        const sw6_bis_Fq XX   = (this->X_).squared();                   // XX  = X1^2
        const sw6_bis_Fq ZZ   = (this->Z_).squared();                   // ZZ  = Z1^2
        const sw6_bis_Fq w    = sw6_bis_G1::coeff_a * ZZ + (XX + XX + XX); // w   = a*ZZ + 3*XX
        const sw6_bis_Fq Y1Z1 = (this->Y_) * (this->Z_);
        const sw6_bis_Fq s    = Y1Z1 + Y1Z1;                            // s   = 2*Y1*Z1
        const sw6_bis_Fq ss   = s.squared();                            // ss  = s^2
        const sw6_bis_Fq sss  = s * ss;                                 // sss = s*ss
        const sw6_bis_Fq R    = (this->Y_) * s;                         // R   = Y1*s
        const sw6_bis_Fq RR   = R.squared();                            // RR  = R^2
        const sw6_bis_Fq B    = ((this->X_)+R).squared()-XX-RR;         // B   = (X1+R)^2 - XX - RR
        const sw6_bis_Fq h    = w.squared() - (B+B);                    // h   = w^2 - 2*B
        const sw6_bis_Fq X3   = h * s;                                  // X3  = h*s
        const sw6_bis_Fq Y3   = w * (B-h)-(RR+RR);                      // Y3  = w*(B-h) - 2*RR
        const sw6_bis_Fq Z3   = sss;                                    // Z3  = sss

        return sw6_bis_G1(X3, Y3, Z3);
    }
}

bool sw6_bis_G1::is_well_formed() const
{
    if (this->is_zero())
    {
        return true;
    }
    else
    {
        /*
          y^2 = x^3 + ax + b

          We are using projective, so equation we need to check is actually

          (y/z)^2 = (x/z)^3 + a (x/z) + b
          z y^2 = x^3  + a z^2 x + b z^3

          z (y^2 - b z^2) = x ( x^2 + a z^2)
        */
        const sw6_bis_Fq X2 = this->X_.squared();
        const sw6_bis_Fq Y2 = this->Y_.squared();
        const sw6_bis_Fq Z2 = this->Z_.squared();

        return (this->Z_ * (Y2 - sw6_bis_G1::coeff_b * Z2) == this->X_ * (X2 + sw6_bis_G1::coeff_a * Z2));
    }
}

sw6_bis_G1 sw6_bis_G1::zero()
{
    return G1_zero;
}

sw6_bis_G1 sw6_bis_G1::one()
{
    return G1_one;
}

sw6_bis_G1 sw6_bis_G1::random_element()
{
    return (scalar_field::random_element().as_bigint()) * G1_one;
}

std::ostream& operator<<(std::ostream &out, const sw6_bis_G1 &g)
{
    sw6_bis_G1 copy(g);
    copy.to_affine_coordinates();

    out << (copy.is_zero() ? 1 : 0) << OUTPUT_SEPARATOR;
#ifdef NO_PT_COMPRESSION
    out << copy.X_ << OUTPUT_SEPARATOR << copy.Y_;
#else
    /* storing LSB of Y */
    out << copy.X_ << OUTPUT_SEPARATOR << (copy.Y_.as_bigint().data[0] & 1);
#endif

    return out;
}

std::istream& operator>>(std::istream &in, sw6_bis_G1 &g)
{
    char is_zero;
    sw6_bis_Fq tX, tY;

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

    // y = +/- sqrt(x^3 + a*x + b)
    if (is_zero)
    {
        sw6_bis_Fq tX2 = tX.squared();
        sw6_bis_Fq tY2 = (tX2 + sw6_bis_G1::coeff_a) * tX + sw6_bis_G1::coeff_b;
        tY = tY2.sqrt();

        if ((tY.as_bigint().data[0] & 1) != Y_lsb)
        {
            tY = -tY;
        }
    }
#endif
    // using projective coordinates
    if (!is_zero)
    {
        g.X_ = tX;
        g.Y_ = tY;
        g.Z_ = sw6_bis_Fq::one();
    }
    else
    {
        g = sw6_bis_G1::zero();
    }

    return in;
}

std::ostream& operator<<(std::ostream& out, const std::vector<sw6_bis_G1> &v)
{
    out << v.size() << "\n";
    for (const sw6_bis_G1& t : v)
    {
        out << t << OUTPUT_NEWLINE;
    }

    return out;
}

std::istream& operator>>(std::istream& in, std::vector<sw6_bis_G1> &v)
{
    v.clear();

    size_t s;
    in >> s;
    consume_newline(in);

    v.reserve(s);

    for (size_t i = 0; i < s; ++i)
    {
        sw6_bis_G1 g;
        in >> g;
        consume_OUTPUT_NEWLINE(in);
        v.emplace_back(g);
    }

    return in;
}

void sw6_bis_G1::batch_to_special_all_non_zeros(std::vector<sw6_bis_G1> &vec)
{
    std::vector<sw6_bis_Fq> Z_vec;
    Z_vec.reserve(vec.size());

    for (auto &el: vec)
    {
        Z_vec.emplace_back(el.Z());
    }
    batch_invert<sw6_bis_Fq>(Z_vec);

    const sw6_bis_Fq one = sw6_bis_Fq::one();

    for (size_t i = 0; i < vec.size(); ++i)
    {
        vec[i] = sw6_bis_G1(vec[i].X() * Z_vec[i], vec[i].Y() * Z_vec[i], one);
    }
}

} // libff
