#pragma once

#include <cmath>
#include "./mathUtils.h"

template <typename T>
class Vec2 {
    public:
         Vec2();
         Vec2(const Vec2<T>& _v);
         Vec2(Vec2<T>&& _rhs) noexcept;
         Vec2(const T& _x, const T& _y);
        ~Vec2() noexcept;

        Vec2<T>& operator=(const Vec2<T>& _v);
        Vec2<T>& operator=(Vec2<T>&& _rhs) noexcept;

        Vec2<T>& operator+=(const Vec2<T>& _v);
        Vec2<T>& operator-=(const Vec2<T>& _v);
        Vec2<T>& operator+=(const T& _val);
        Vec2<T>& operator-=(const T& _val);
        Vec2<T>& operator*=(const T& _val);
        Vec2<T>& operator/=(const T& _val);

        Vec2<T> operator+(const Vec2<T>& _v);
        Vec2<T> operator-(const Vec2<T>& _v);
        Vec2<T> operator+(const T& _val);
        Vec2<T> operator-(const T& _val);
        Vec2<T> operator*(const T& _val);
        Vec2<T> operator/(const T& _val);

        const T cross(const Vec2<T>& _v) const;
        const T dot  (const Vec2<T>& _v) const;

        const double magnitude  () const;
        const float  magnitude_f() const;

        const double distance  (const Vec2<T>& _v) const;
        const float  distance_f(const Vec2<T>& _v) const;

        Vec2<T>      normalize  () const;
        Vec2<float>  normalize_f() const;
        Vec2<double> normalize_d() const;

        const float angle(const Vec2<T>& _v) const;

        void x(const T& _x);
        void y(const T& _y);

        constexpr T x() const;
        constexpr T y() const;

        const T* ptr() const;

        static const T cross(const Vec2<T>& _v1, const Vec2<T>& _v2);
        static const T dot  (const Vec2<T>& _v1, const Vec2<T>& _v2);

        static const double magnitude  (const Vec2<T>& _v);
        static const float  magnitude_f(const Vec2<T>& _v);

        static const double distance  (const Vec2<T>& _v1, const Vec2<T>& _v2);
        static const float  distance_f(const Vec2<T>& _v1, const Vec2<T>& _v2);

        static Vec2<T>      normalize  (const Vec2<T>& _v);
        static Vec2<float>  normalize_f(const Vec2<T>& _v);
        static Vec2<double> normalize_d(const Vec2<T>& _v);

        static const float angle(const Vec2<T>& _v1, const Vec2<T>& _v2);

    private:
        T m_val[2];
};

template <typename T> Vec2<T>::Vec2()
    : m_val{} { }
template <typename T> Vec2<T>::Vec2(const Vec2<T>& _v)
    : m_val{ _v.m_val[VEC::X], _v.m_val[VEC::Y] } { }
template <typename T> Vec2<T>::Vec2(Vec2<T>&& _rhs) noexcept
    : m_val{ _rhs.m_val[VEC::X], _rhs.m_val[VEC::Y] } {
    _rhs.m_val[VEC::X] = 0;
    _rhs.m_val[VEC::Y] = 0;
}
template <typename T> Vec2<T>::Vec2(const T& _x, const T& _y)
    : m_val{ _x, _y } { }
template <typename T> Vec2<T>::~Vec2() noexcept { }

template <typename T> Vec2<T>& Vec2<T>::operator=(const Vec2<T>& _v) {
    m_val[VEC::X] = _v.m_val[VEC::X];
    m_val[VEC::Y] = _v.m_val[VEC::Y];

    return *this;
}
template <typename T> Vec2<T>& Vec2<T>::operator=(Vec2<T>&& _rhs) noexcept {
    m_val[VEC::X] = _rhs.m_val[VEC::X];
    m_val[VEC::Y] = _rhs.m_val[VEC::Y];

    _rhs.m_val[VEC::X] = 0;
    _rhs.m_val[VEC::Y] = 0;

    return *this;
}

template <typename T> Vec2<T>& Vec2<T>::operator+=(const Vec2<T>& _v) {
    m_val[VEC::X] += _v.m_val[VEC::X];
    m_val[VEC::Y] += _v.m_val[VEC::Y];

    return *this;
}
template <typename T> Vec2<T>& Vec2<T>::operator-=(const Vec2<T>& _v) {
    m_val[VEC::X] -= _v.m_val[VEC::X];
    m_val[VEC::Y] -= _v.m_val[VEC::Y];

    return *this;
}
template <typename T> Vec2<T>& Vec2<T>::operator+=(const T& _val) {
    m_val[VEC::X] += _val;
    m_val[VEC::Y] += _val;

    return *this;
}
template <typename T> Vec2<T>& Vec2<T>::operator-=(const T& _val) {
    m_val[VEC::X] -= _val;
    m_val[VEC::Y] -= _val;

    return *this;
}
template <typename T> Vec2<T>& Vec2<T>::operator*=(const T& _val) {
    m_val[VEC::X] *= _val;
    m_val[VEC::Y] *= _val;

    return *this;
}
template <typename T> Vec2<T>& Vec2<T>::operator/=(const T& _val) {
    m_val[VEC::X] /= _val;
    m_val[VEC::Y] /= _val;

    return *this;
}

template <typename T> Vec2<T> Vec2<T>::operator+(const Vec2<T>& _v) {
    return Vec2<T>(
        m_val[VEC::X] + _v.m_val[VEC::X],
        m_val[VEC::Y] + _v.m_val[VEC::Y]
    );
}
template <typename T> Vec2<T> Vec2<T>::operator-(const Vec2<T>& _v) {
    return Vec2<T>(
        m_val[VEC::X] - _v.m_val[VEC::X],
        m_val[VEC::Y] - _v.m_val[VEC::Y]
    );
}
template <typename T> Vec2<T> Vec2<T>::operator+(const T& _val)     {
    return Vec2<T>(
        m_val[VEC::X] + _val,
        m_val[VEC::Y] + _val
    );
}
template <typename T> Vec2<T> Vec2<T>::operator-(const T& _val)     {
    return Vec2<T>(
        m_val[VEC::X] - _val,
        m_val[VEC::Y] - _val
    );
}
template <typename T> Vec2<T> Vec2<T>::operator*(const T& _val)     {
    return Vec2<T>(
        m_val[VEC::X] * _val,
        m_val[VEC::Y] * _val
    );
}
template <typename T> Vec2<T> Vec2<T>::operator/(const T& _val)     {
    return Vec2<T>(
        m_val[VEC::X] / _val,
        m_val[VEC::Y] / _val
    );
}

template <typename T> const T Vec2<T>::cross(const Vec2<T>& _v) const {
    return (
        m_val[VEC::X] * _v.m_val[VEC::Y] -
        m_val[VEC::Y] * _v.m_val[VEC::X]
    );
}
template <typename T> const T Vec2<T>::dot  (const Vec2<T>& _v) const {
    return (
        m_val[VEC::X] * _v.m_val[VEC::X] +
        m_val[VEC::Y] * _v.m_val[VEC::Y]
    );
}

template <typename T> const double Vec2<T>::magnitude  () const {
    return sqrt(
        MATH::SQUARE(m_val[VEC::X]) +
        MATH::SQUARE(m_val[VEC::Y])
    );
}
template <typename T> const float  Vec2<T>::magnitude_f() const {
    return sqrtf(
        MATH::SQUARE(m_val[VEC::X]) +
        MATH::SQUARE(m_val[VEC::Y])
    );
}

template <typename T> const double Vec2<T>::distance  (const Vec2<T>& _v) const { 
   return sqrt(
       MATH::SQUARE(m_val[VEC::X] - _v.m_val[VEC::X]) +
       MATH::SQUARE(m_val[VEC::Y] - _v.m_val[VEC::Y])
   );
}
template <typename T> const float  Vec2<T>::distance_f(const Vec2<T>& _v) const {
    return sqrtf(
        MATH::SQUARE(m_val[VEC::X] - _v.m_val[VEC::X]) +
        MATH::SQUARE(m_val[VEC::Y] - _v.m_val[VEC::Y])
    );
}

template <typename T> Vec2<T>      Vec2<T>::normalize  () const {
    const float mag = magnitude_f();

    return Vec2<T>((T)(m_val[VEC::X] / mag), (T)(m_val[VEC::Y] / mag));
}
template <typename T> Vec2<float>  Vec2<T>::normalize_f() const {
    const float mag = magnitude_f();

    return (Vec2<float>(m_val[VEC::X], m_val[VEC::Y]) /= mag);
}
template <typename T> Vec2<double> Vec2<T>::normalize_d() const {
    const double mag = magnitude();

    return (Vec2<double>(m_val[VEC::X], m_val[VEC::Y]) /= mag);
}

template <typename T> const float Vec2<T>::angle(const Vec2<T>& _v) const { return MATH::RAD_TO_DEG(acosf((dot(_v) / (magnitude_f() * _v.magnitude_f())))); }

template <typename T> void Vec2<T>::x(const T& _x) { m_val[VEC::X] = _x; }
template <typename T> void Vec2<T>::y(const T& _y) { m_val[VEC::Y] = _y; }

template <typename T> constexpr T Vec2<T>::x() const { return m_val[VEC::X]; }
template <typename T> constexpr T Vec2<T>::y() const { return m_val[VEC::Y]; }

template <typename T> const T* Vec2<T>::ptr() const { return m_val; }

template <typename T> const T Vec2<T>::cross(const Vec2<T>& _v1, const Vec2<T>& _v2) { return _v1.cross(_v2); }
template <typename T> const T Vec2<T>::dot  (const Vec2<T>& _v1, const Vec2<T>& _v2) { return _v1.dot(_v2);   }

template <typename T> const double Vec2<T>::magnitude  (const Vec2<T>& _v) { return _v.magnitude();   }
template <typename T> const float  Vec2<T>::magnitude_f(const Vec2<T>& _v) { return _v.magnitude_f(); }

template <typename T> const double Vec2<T>::distance  (const Vec2<T>& _v1, const Vec2<T>& _v2) { return _v1.distance(_v2);   }
template <typename T> const float  Vec2<T>::distance_f(const Vec2<T>& _v1, const Vec2<T>& _v2) { return _v1.distance_f(_v2); }

template <typename T> Vec2<T>      Vec2<T>::normalize  (const Vec2<T>& _v) { return _v.normalize();   }
template <typename T> Vec2<float>  Vec2<T>::normalize_f(const Vec2<T>& _v) { return _v.normalize_f(); }
template <typename T> Vec2<double> Vec2<T>::normalize_d(const Vec2<T>& _v) { return _v.normalize_d(); }

template <typename T> const float Vec2<T>::angle(const Vec2<T>& _v1, const Vec2<T>& _v2) { return _v1.angle(_v2); }