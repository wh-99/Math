#pragma once

#include <cmath>
#include "./mathUtils.h"

template <typename T>
class Vec3 {
    public:
         Vec3();
         Vec3(const Vec3<T>& _v);
         Vec3(Vec3<T>&& _rhs) noexcept;
         Vec3(const T& _x, const T& _y, const T& _z);
        ~Vec3() noexcept;

        Vec3<T>& operator=(const Vec3<T>& _v);
        Vec3<T>& operator=(Vec3<T>&& _rhs) noexcept;
        
        Vec3<T>& operator+=(const Vec3<T>& _v);
        Vec3<T>& operator-=(const Vec3<T>& _v);
        Vec3<T>& operator+=(const T& _val);
        Vec3<T>& operator-=(const T& _val);
        Vec3<T>& operator*=(const T& _val);
        Vec3<T>& operator/=(const T& _val);

        Vec3<T> operator+(const Vec3<T>& _v);
        Vec3<T> operator-(const Vec3<T>& _v);
        Vec3<T> operator+(const T& _val);
        Vec3<T> operator-(const T& _val);
        Vec3<T> operator*(const T& _val);
        Vec3<T> operator/(const T& _val);

        Vec3<T> cross(const Vec3<T>& _v) const;
        const T dot  (const Vec3<T>& _v) const;

        const double magnitude  () const;
        const float  magnitude_f() const;

        const double distance  (const Vec3<T>& _v) const;
        const float  distance_f(const Vec3<T>& _v) const;

        Vec3<T>      normalize  () const;
        Vec3<float>  normalize_f() const;
        Vec3<double> normalize_d() const;

        const float angle(const Vec3<T>& _v) const;

        void x(const T& _x);
        void y(const T& _y);
        void z(const T& _z);

        constexpr T x() const;
        constexpr T y() const;
        constexpr T z() const;

        const T* ptr() const;

        static Vec3<T> cross(const Vec3<T>& _v1, const Vec3<T>& _v2);
        static const T dot  (const Vec3<T>& _v1, const Vec3<T>& _v2);

        static const double magnitude  (const Vec3<T>& _v);
        static const float  magnitude_f(const Vec3<T>& _v);

        static const double distance  (const Vec3<T>& _v1, const Vec3<T>& _v2);
        static const float  distance_f(const Vec3<T>& _v1, const Vec3<T>& _v2);

        static Vec3<T>      normalize  (const Vec3<T>& _v);
        static Vec3<float>  normalize_f(const Vec3<T>& _v);
        static Vec3<double> normalize_d(const Vec3<T>& _v);

        static const float angle(const Vec3<T>& _v1, const Vec3<T>& _v2);

    private:
        T m_val[3];
};

template <typename T> Vec3<T>::Vec3()
    : m_val{} { }
template <typename T> Vec3<T>::Vec3(const Vec3<T>& _v)
    : m_val{ _v.m_val[VEC::X], _v.m_val[VEC::Y], _v.m_val[VEC::Z] } { }
template <typename T> Vec3<T>::Vec3(Vec3<T>&& _rhs) noexcept
    : m_val{ _rhs.m_val[VEC::X], _rhs.m_val[VEC::Y], _rhs.m_val[VEC::Z] } {
    _rhs.m_val[VEC::X] = 0;
    _rhs.m_val[VEC::Y] = 0;
    _rhs.m_val[VEC::Z] = 0;
}
template <typename T> Vec3<T>::Vec3(const T& _x, const T& _y, const T& _z)
    : m_val{ _x, _y, _z } { }
template <typename T> Vec3<T>::~Vec3() noexcept { }

template <typename T> Vec3<T>& Vec3<T>::operator=(const Vec3<T>& _v) {
    m_val[VEC::X] = _v.m_val[VEC::X];
    m_val[VEC::Y] = _v.m_val[VEC::Y];
    m_val[VEC::Z] = _v.m_val[VEC::Z];

    return *this;
}
template <typename T> Vec3<T>& Vec3<T>::operator=(Vec3<T>&& _rhs) noexcept {
    m_val[VEC::X] = _rhs.m_val[VEC::X];
    m_val[VEC::Y] = _rhs.m_val[VEC::Y];
    m_val[VEC::Z] = _rhs.m_val[VEC::Z];

    _rhs.m_val[VEC::X] = 0;
    _rhs.m_val[VEC::Y] = 0;
    _rhs.m_val[VEC::Z] = 0;

    return *this;
}

template <typename T> Vec3<T>& Vec3<T>::operator+=(const Vec3<T>& _v) {
    m_val[VEC::X] += _v.m_val[VEC::X];
    m_val[VEC::Y] += _v.m_val[VEC::Y];
    m_val[VEC::Z] += _v.m_val[VEC::Z];

    return *this;
}
template <typename T> Vec3<T>& Vec3<T>::operator-=(const Vec3<T>& _v) {
    m_val[VEC::X] -= _v.m_val[VEC::X];
    m_val[VEC::Y] -= _v.m_val[VEC::Y];
    m_val[VEC::Z] -= _v.m_val[VEC::Z];

    return *this;
}
template <typename T> Vec3<T>& Vec3<T>::operator+=(const T& _val) {
    m_val[VEC::X] += _val;
    m_val[VEC::Y] += _val;
    m_val[VEC::Z] += _val;

    return *this;
}
template <typename T> Vec3<T>& Vec3<T>::operator-=(const T& _val) {
    m_val[VEC::X] -= _val;
    m_val[VEC::Y] -= _val;
    m_val[VEC::Z] -= _val;

    return *this;
}
template <typename T> Vec3<T>& Vec3<T>::operator*=(const T& _val) {
    m_val[VEC::X] *= _val;
    m_val[VEC::Y] *= _val;
    m_val[VEC::Z] *= _val;

    return *this;
}
template <typename T> Vec3<T>& Vec3<T>::operator/=(const T& _val) {
    m_val[VEC::X] /= _val;
    m_val[VEC::Y] /= _val;
    m_val[VEC::Z] /= _val;

    return *this;
}

template <typename T> Vec3<T> Vec3<T>::operator+(const Vec3<T>& _v) {
    return Vec3<T>(
        m_val[VEC::X] + _v.m_val[VEC::X],
        m_val[VEC::Y] + _v.m_val[VEC::Y],
        m_val[VEC::Z] + _v.m_val[VEC::Z]
    );
}
template <typename T> Vec3<T> Vec3<T>::operator-(const Vec3<T>& _v) {
    return Vec3<T>(
        m_val[VEC::X] - _v.m_val[VEC::X],
        m_val[VEC::Y] - _v.m_val[VEC::Y],
        m_val[VEC::Z] - _v.m_val[VEC::Z]
    );
}
template <typename T> Vec3<T> Vec3<T>::operator+(const T& _val)     {
    return Vec3<T>(
        m_val[VEC::X] + _val,
        m_val[VEC::Y] + _val,
        m_val[VEC::Z] + _val
    );
}
template <typename T> Vec3<T> Vec3<T>::operator-(const T& _val)     {
    return Vec3<T>(
        m_val[VEC::X] - _val,
        m_val[VEC::Y] - _val,
        m_val[VEC::Z] - _val
    );
}
template <typename T> Vec3<T> Vec3<T>::operator*(const T& _val)     {
    return Vec3<T>(
        m_val[VEC::X] * _val,
        m_val[VEC::Y] * _val,
        m_val[VEC::Z] * _val
    );
}
template <typename T> Vec3<T> Vec3<T>::operator/(const T& _val)     {
    return Vec3<T>(
        m_val[VEC::X] / _val,
        m_val[VEC::Y] / _val,
        m_val[VEC::Z] / _val
    );
}

template <typename T> Vec3<T> Vec3<T>::cross(const Vec3<T>& _v) const {
    return Vec3<T>(
        m_val[VEC::Y] * _v.m_val[VEC::Z] - m_val[VEC::Z] * _v.m_val[VEC::Y],
        m_val[VEC::Z] * _v.m_val[VEC::X] - m_val[VEC::X] * _v.m_val[VEC::Z],
        m_val[VEC::X] * _v.m_val[VEC::Y] - m_val[VEC::Y] * _v.m_val[VEC::X]
    );
}
template <typename T> const T Vec3<T>::dot  (const Vec3<T>& _v) const {
    return (
        m_val[VEC::X] * _v.m_val[VEC::X] +
        m_val[VEC::Y] * _v.m_val[VEC::Y] +
        m_val[VEC::Z] * _v.m_val[VEC::Z]
    );
}

template <typename T> const double Vec3<T>::magnitude  () const {
    return sqrt(
        MATH::SQUARE(m_val[VEC::X]) +
        MATH::SQUARE(m_val[VEC::Y]) +
        MATH::SQUARE(m_val[VEC::Z])
    );
}
template <typename T> const float  Vec3<T>::magnitude_f() const {
    return sqrtf(
        MATH::SQUARE(m_val[VEC::X]) +
        MATH::SQUARE(m_val[VEC::Y]) +
        MATH::SQUARE(m_val[VEC::Z])
    );
}

template <typename T> const double Vec3<T>::distance  (const Vec3<T>& _v) const {
    return sqrt(
        MATH::SQUARE(m_val[VEC::X] - _v.m_val[VEC::X]) +
        MATH::SQUARE(m_val[VEC::Y] - _v.m_val[VEC::Y]) +
        MATH::SQUARE(m_val[VEC::Z] - _v.m_val[VEC::Z])
    );
}
template <typename T> const float  Vec3<T>::distance_f(const Vec3<T>& _v) const {
    return sqrtf(
        MATH::SQUARE(m_val[VEC::X] - _v.m_val[VEC::X]) +
        MATH::SQUARE(m_val[VEC::Y] - _v.m_val[VEC::Y]) +
        MATH::SQUARE(m_val[VEC::Z] - _v.m_val[VEC::Z])
    );
}

template <typename T> Vec3<T>      Vec3<T>::normalize  () const {
    float mag = magnitude_f();

    Vec3<T>((T)(m_val[VEC::X] / mag), (T)(m_val[VEC::Y] / mag), (T)(m_val[VEC::Z] / mag));
}
template <typename T> Vec3<float>  Vec3<T>::normalize_f() const {
    float mag = magnitude_f();

    return (Vec3<float>(m_val[VEC::X], m_val[VEC::Y], m_val[VEC::Z]) /= mag);
}
template <typename T> Vec3<double> Vec3<T>::normalize_d() const {
    double mag = magnitude();

    return (Vec3<double>(m_val[VEC::X], m_val[VEC::Y], m_val[VEC::Z]) /= mag);
}

template <typename T> const float Vec3<T>::angle(const Vec3<T>& _v) const { return MATH::RAD_TO_DEG(acosf((dot(_v) / (magnitude_f() * _v.magnitude_f())))); }

template <typename T> void Vec3<T>::x(const T& _x) { m_val[VEC::X] = _x; }
template <typename T> void Vec3<T>::y(const T& _y) { m_val[VEC::Y] = _y; }
template <typename T> void Vec3<T>::z(const T& _z) { m_val[VEC::Z] = _z; }

template <typename T> constexpr T Vec3<T>::x() const { return m_val[VEC::X]; }
template <typename T> constexpr T Vec3<T>::y() const { return m_val[VEC::Y]; }
template <typename T> constexpr T Vec3<T>::z() const { return m_val[VEC::Z]; }

template <typename T> const T* Vec3<T>::ptr() const { return m_val; }

template <typename T> Vec3<T> Vec3<T>::cross(const Vec3<T>& _v1, const Vec3<T>& _v2) { return _v1.cross(_v2); }
template <typename T> const T Vec3<T>::dot  (const Vec3<T>& _v1, const Vec3<T>& _v2) { return _v1.dot(_v2);   }

template <typename T> const double Vec3<T>::magnitude  (const Vec3<T>& _v) { return _v.magnitude();   }
template <typename T> const float  Vec3<T>::magnitude_f(const Vec3<T>& _v) { return _v.magnitude_f(); }

template <typename T> const double Vec3<T>::distance  (const Vec3<T>& _v1, const Vec3<T>& _v2) { return _v1.distance(_v2);   }
template <typename T> const float  Vec3<T>::distance_f(const Vec3<T>& _v1, const Vec3<T>& _v2) { return _v1.distance_f(_v2); }

template <typename T> Vec3<T>      Vec3<T>::normalize  (const Vec3<T>& _v) { return _v.normalize();   }
template <typename T> Vec3<float>  Vec3<T>::normalize_f(const Vec3<T>& _v) { return _v.normalize_f(); }
template <typename T> Vec3<double> Vec3<T>::normalize_d(const Vec3<T>& _v) { return _v.normalize_d(); }

template <typename T> const float Vec3<T>::angle(const Vec3<T>& _v1, const Vec3<T>& _v2) { return _v1.angle(_v2); }