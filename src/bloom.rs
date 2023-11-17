use ahash::RandomState;
use bit_vec::BitVec;
use core::hash::Hash;

pub struct BloomFilter {
    size: usize,
    n_hashes: usize,
    bv: BitVec,
    hash_builders: (RandomState, RandomState),
}

impl BloomFilter {
    const BLOCK_SIZE: usize = 1 << 12;
    const BLOCK_MASK: usize = Self::BLOCK_SIZE - 1;
    const BLOCK_PREFIX: usize = !Self::BLOCK_MASK;

    pub fn new_with_seed(size: usize, n_hashes: usize, seed: u64) -> Self {
        let size = size.saturating_add(Self::BLOCK_SIZE - 1) / Self::BLOCK_SIZE * Self::BLOCK_SIZE;
        Self {
            size,
            n_hashes,
            bv: BitVec::from_elem(size, false),
            hash_builders: (
                RandomState::with_seeds(seed, seed + 1, seed + 2, seed + 3),
                RandomState::with_seeds(seed + 4, seed + 5, seed + 6, seed + 7),
            ),
        }
    }

    pub fn new(size: usize, n_hashes: usize) -> Self {
        let seed = (size + n_hashes) as u64;
        Self::new_with_seed(size, n_hashes, seed)
    }

    fn hashes<T: Hash>(&self, x: T) -> (u64, u64) {
        (
            self.hash_builders.0.hash_one(&x),
            self.hash_builders.1.hash_one(&x),
        )
    }

    fn indices<T: Hash>(&self, x: T) -> Vec<usize> {
        let mut res = Vec::with_capacity(self.n_hashes);
        let (h0, h1) = self.hashes(x);
        let u = h0 as usize % self.size;
        let v = h1 as usize;
        let block_addr = u & Self::BLOCK_PREFIX;
        let mut local_addr = u;
        res.push(u);
        (1..self.n_hashes).for_each(|_| {
            local_addr = (local_addr + v) & Self::BLOCK_MASK;
            res.push(block_addr | local_addr);
        });
        res
    }

    pub fn contains<T: Hash>(&self, x: T) -> bool {
        self.indices(x)
            .iter()
            .all(|&i| self.bv.get(i).unwrap_or(false))
    }

    pub fn insert<T: Hash>(&mut self, x: T) {
        self.indices(x).iter().for_each(|&i| self.bv.set(i, true));
    }

    pub fn insert_if_missing<T: Hash>(&mut self, x: T) -> bool {
        let mut missing = false;
        for i in self.indices(x) {
            if !self.bv.get(i).unwrap_or(false) {
                missing = true;
                self.bv.set(i, true);
            }
        }
        missing
    }

    pub fn clear(&mut self){
        self.bv.clear();
    }
}
