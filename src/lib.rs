//! This module provides the Gf256 type which is used to represent
//! elements of a finite field wich 256 elements.
//! Based on 
//! https://github.com/sellibitze/secretshare/blob/
//!     master/src/gf256.rs

use std::num::Wrapping;
use std::ops::{ Add, Sub, Mul, Div };
use std::sync::{ Once, ONCE_INIT };

const POLY: u8 = 0x1D; // represents x^8 + x^4 + x^3 + x^2 + 1

/// replicates the least significant bit to every other bit
#[inline]
fn mask(bit: u8) -> u8 {
    (Wrapping(0u8) - Wrapping(bit & 1)).0
}

/// multiplies a polynomial with x and returns the residual
/// of the polynomial division with POLY as divisor
#[inline]
fn xtimes(poly: u8) -> u8 {
	(poly << 1) ^ (mask(poly >> 7) & POLY)
}

/// Tables used for multiplication and division
struct Tables {
	exp: [u8; 256],
	log: [u8; 256],
	inv: [u8; 256]
}

static INIT: Once = ONCE_INIT;
static mut TABLES: Tables = Tables {
	exp: [0; 256],
	log: [0; 256],
	inv: [0; 256]
};

fn get_tables() -> &'static Tables {
	INIT.call_once(|| {
		// mutable access is fine because of synchronization via INIT
		let tabs = unsafe { &mut TABLES };
		let mut tmp = 1;
		for power in 0..255usize {
			tabs.exp[power] = tmp;
			tabs.log[tmp as usize] = power as u8;
			tmp = xtimes(tmp);
		}
		tabs.exp[255] = 1;
		for x in 1..256usize {
			let l = tabs.log[x];
			let nl = if l == 0 { 0 } else { 255 - l };
			let i = tabs.exp[nl as usize];
			tabs.inv[x] = i;
		}
	});
	// We're guaranteed to have TABLES initialized by now
	return unsafe { &TABLES };
}

/// Type for elements of a finite field with 256 elements
#[derive(Copy,Clone,PartialEq,Eq,Debug)]
pub struct Gf256 {
	pub poly: u8
}

impl Gf256 {
	/// returns the additive neutral element of the field
	#[inline]
	pub fn zero() -> Gf256 {
		Gf256 { poly: 0 }
	}
	/// returns the multiplicative neutral element of the field
	#[inline]
	pub fn one() -> Gf256 {
		Gf256 { poly: 1 }
	}
	#[inline]
	pub fn from_byte(b: u8) -> Gf256 {
		Gf256 { poly: b }
	}
	#[inline]
	pub fn to_byte(&self) -> u8 {
		self.poly
	}

    /// Find the logrithm of self in base x.
    /// In other words, find a value i such that x ^ i == self
	pub fn log(&self) -> Option<u8> {
		if self.poly == 0 {
			None
		} else {
			let tabs = get_tables();
			Some(tabs.log[self.poly as usize])
		}
	}

    /// Calculate x ^ power
	pub fn xexp(power: u8) -> Gf256 {
		let tabs = get_tables();
		Gf256 { poly: tabs.exp[power as usize] }
	}

    /// Calculate self ^ power
	pub fn exp(&self, power: u8) -> Gf256 {
        match self.log() {
            None => Self::zero(),
            Some(i) => { 
                // Current value is x^i
                Self::xexp((((i as u16) * (power as u16)) % 255) as u8)
            },
        }
	}

    /// Find the inverse of self: A number y such that self * y == 1
	pub fn inv(&self) -> Option<Gf256> {
		self.log().map(|l| Gf256::xexp(255 - l))
	}
}

impl Add<Gf256> for Gf256 {
	type Output = Gf256;
	#[inline]
	fn add(self, rhs: Gf256) -> Gf256 {
		Gf256::from_byte(self.poly ^ rhs.poly)
	}
}

impl Sub<Gf256> for Gf256 {
	type Output = Gf256;
	#[inline]
	fn sub(self, rhs: Gf256) -> Gf256 {
		Gf256::from_byte(self.poly ^ rhs.poly)
	}
}

impl Mul<Gf256> for Gf256 {
	type Output = Gf256;
	fn mul(self, rhs: Gf256) -> Gf256 {
		if let (Some(l1), Some(l2)) = (self.log(), rhs.log()) {
			let tmp = ((l1 as u16) + (l2 as u16)) % 255;
			Gf256::xexp(tmp as u8)
		} else {
			Gf256 { poly: 0 }
		}
	}
}

impl Div<Gf256> for Gf256 {
	type Output = Gf256;
	fn div(self, rhs: Gf256) -> Gf256 {
		let l2 = rhs.log().expect("Division by zero");
		if let Some(l1) = self.log() {
			let tmp = ((l1 as u16) + 255 - (l2 as u16)) % 255;
			Gf256::xexp(tmp as u8)
		} else {
			Gf256 { poly: 0 }
		}
	}
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_sub() {
        let a = Gf256::from_byte(0x5);
        let b = Gf256::from_byte(0xab);

        let c = a + b;

        let b_tag = c - a;
        assert_eq!(b,b_tag);

        let a_tag = c - b;
        assert_eq!(a,a_tag);
    }

    #[test]
    fn test_mul_div() {
        let a = Gf256::from_byte(0x5);
        let b = Gf256::from_byte(0xab);

        let c = a * b;

        let b_tag = c / a;
        assert_eq!(b,b_tag);

        let a_tag = c / b;
        assert_eq!(a,a_tag);
    }

    #[test]
    #[should_panic]
    fn test_div_by_0() {
        let a = Gf256::from_byte(0x5);
        let b = Gf256::from_byte(0);

        // Attempt to divide by zero. This should panic:
        let _ = a / b;
    }

    #[test]
    fn test_distributive() {
        let a = Gf256::from_byte(0x8f);
        let b = Gf256::from_byte(0x15);
        let c = Gf256::from_byte(0xa2);

        assert_eq!(a * (b + c), a*b + a*c);
    }

    #[test]
    fn test_commutative() {
        let a = Gf256::from_byte(0x8f);
        let b = Gf256::from_byte(0x15);

        assert_eq!(a + b, b + a);
        assert_eq!(a * b, b * a);
    }

    #[test]
    fn test_associative() {
        let a = Gf256::from_byte(0x8f);
        let b = Gf256::from_byte(0x15);
        let c = Gf256::from_byte(0xa2);

        assert_eq!((a + b) + c, a + (b + c));
        assert_eq!((a * b) * c, a * (b * c));
    }

    #[test]
    fn test_inverse() {
        let a = Gf256::from_byte(0x8f);
        let b = Gf256::from_byte(0x15);
        let one = Gf256::one();

        assert_eq!(a * a.inv().unwrap(), one);
        assert_eq!(a.inv().unwrap() * a, one);
        assert_eq!(b * b.inv().unwrap(), one);
        assert_eq!(b.inv().unwrap() * b, one);
    }

    #[test]
    fn test_inverse_zero() {
        let zero = Gf256::zero();

        let res = zero.inv();
        assert_eq!(res, None);
    }

    #[test]
    fn test_exp() {
        let a = Gf256::from_byte(0x8f);

        assert_eq!(a.exp(0), Gf256::one());
        assert_eq!(a.exp(1), a);
        assert_eq!(a.exp(2), a*a);
        assert_eq!(a.exp(3), a*a*a);

        // Fermat's little theorem:
        assert_eq!(a.exp(255), Gf256::one());
    }
}

