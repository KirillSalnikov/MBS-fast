#ifndef NUMBERLIKEARRAY_H
#define NUMBERLIKEARRAY_H

// Make sure we have NULL.
#ifndef NULL
#define NULL 0
#endif

/* A NumberlikeArray<Blk> object holds a variable-sized array of Blk with a
 * length and a capacity and provides basic memory management features.
 * BigUnsigned and BigUnsignedInABase both subclass it.
 *
 * NumberlikeArray provides no information hiding.  Subclasses should use
 * nonpublic inheritance and manually expose members as desired using
 * declarations like this:
 *
 * public:
 *     NumberlikeArray< the-type-argument >::getLength;
 */
template <class Blk>
class NumberlikeArray {
public:

	// Type for the index of a block in the array
	typedef unsigned int Index;
	// The number of bits in a block, defined below.
	static const unsigned int N;

	// The current allocated capacity of this NumberlikeArray (in blocks)
	Index cap;
	// The actual length of the value stored in this NumberlikeArray (in blocks)
	Index len;
	// Active storage (can be NULL if len == 0)
	Blk *blk;
private:
	// A short facet history fits in two words on typical convex particles.
	// Larger numbers retain arbitrary precision and the reported capacity.
	Blk inlineBlocks[2];

	Blk *storage(Index count) {
		if (count == 0)
			return NULL;
		if (count <= 2)
			return inlineBlocks;
		return new Blk[count];
	}

	void release(Blk *pointer) {
		if (pointer != inlineBlocks)
			delete [] pointer;
	}
public:

	// Constructs a ``zero'' NumberlikeArray with the given capacity.
	NumberlikeArray(Index c) : cap(c), len(0) { 
		blk = storage(cap);
	}

	// Start with zero active blocks.
	NumberlikeArray() : cap(0), len(0) {
		blk = NULL;
	}

	// Only overflow storage is heap allocated.
	~NumberlikeArray() {
		release(blk);
	}

	/* Ensures that the array has at least the requested capacity; may
	 * destroy the contents. */
	void allocate(Index c);

	/* Ensures that the array has at least the requested capacity; does not
	 * destroy the contents. */
	void allocateAndCopy(Index c);

	// Copy constructor
	NumberlikeArray(const NumberlikeArray<Blk> &x);

	// Assignment operator
	void operator=(const NumberlikeArray<Blk> &x);

	// Constructor that copies from a given array of blocks
	NumberlikeArray(const Blk *b, Index blen);

	// ACCESSORS
	Index getCapacity()     const { return cap;      }
	Index getLength()       const { return len;      }
	Blk   getBlock(Index i) const { return blk[i];   }
	bool  isEmpty()         const { return len == 0; }

	/* Equality comparison: checks if both objects have the same length and
	 * equal (==) array elements to that length.  Subclasses may wish to
	 * override. */
	bool operator ==(const NumberlikeArray<Blk> &x) const;

	bool operator !=(const NumberlikeArray<Blk> &x) const {
		return !operator ==(x);
	}
};

/* BEGIN TEMPLATE DEFINITIONS.  They are present here so that source files that
 * include this header file can generate the necessary real definitions. */

template <class Blk>
const unsigned int NumberlikeArray<Blk>::N = 8 * sizeof(Blk);

template <class Blk>
void NumberlikeArray<Blk>::allocate(Index c) {
	// If the requested capacity is more than the current capacity...
	if (c > cap) {
		// Acquire the new storage before releasing the old value.
		Blk *next = storage(c);
		release(blk);
		blk = next;
		cap = c;
	}
}

template <class Blk>
void NumberlikeArray<Blk>::allocateAndCopy(Index c) {
	// If the requested capacity is more than the current capacity...
	if (c > cap) {
		Blk *oldBlk = blk;
		// Allocate the new number array
		blk = storage(c);
		cap = c;
		// Copy number blocks
		Index i;
		for (i = 0; i < len; i++)
			blk[i] = oldBlk[i];
		// Delete the old array
		release(oldBlk);
	}
}

template <class Blk>
NumberlikeArray<Blk>::NumberlikeArray(const NumberlikeArray<Blk> &x)
		: len(x.len) {
	// Create array
	cap = len;
	blk = storage(cap);
	// Copy blocks
	Index i;
	for (i = 0; i < len; i++)
		blk[i] = x.blk[i];
}

template <class Blk>
void NumberlikeArray<Blk>::operator=(const NumberlikeArray<Blk> &x) {
	/* Calls like a = a have no effect; catch them before the aliasing
	 * causes a problem */
	if (this == &x)
		return;
	allocate(x.len);
	len = x.len;
	// Copy number blocks
	Index i;
	for (i = 0; i < len; i++)
		blk[i] = x.blk[i];
}

template <class Blk>
NumberlikeArray<Blk>::NumberlikeArray(const Blk *b, Index blen)
		: cap(blen), len(blen) {
	// Create array
	blk = storage(cap);
	// Copy blocks
	Index i;
	for (i = 0; i < len; i++)
		blk[i] = b[i];
}

template <class Blk>
bool NumberlikeArray<Blk>::operator ==(const NumberlikeArray<Blk> &x) const {
	if (len != x.len)
		// Definitely unequal.
		return false;
	else {
		// Compare corresponding blocks one by one.
		Index i;
		for (i = 0; i < len; i++)
			if (blk[i] != x.blk[i])
				return false;
		// No blocks differed, so the objects are equal.
		return true;
	}
}

#endif
