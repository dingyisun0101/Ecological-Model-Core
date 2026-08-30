use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};

use ecological_state_toolkit::trajectory::{
    AbundanceView, EquilibriumEvidence, TerminalPolicy, TrajectoryObservation,
    TrajectoryObservationPolicy, TrajectoryObserver,
};

struct CountingAllocator;
static ALLOCATIONS: AtomicUsize = AtomicUsize::new(0);

unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        // SAFETY: the valid layout is forwarded unchanged to the system allocator.
        unsafe { System.alloc(layout) }
    }

    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        // SAFETY: the pointer and layout received from GlobalAlloc are forwarded unchanged.
        unsafe { System.dealloc(pointer, layout) }
    }
}

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

#[test]
fn terminal_observation_allocates_nothing_after_arena_warmup() {
    let mut observer = TrajectoryObserver::from_policy(TrajectoryObservationPolicy::TerminalOnly(
        TerminalPolicy {
            sample_interval_iterations: 1,
            trailing_window_samples: 3,
        },
    ))
    .unwrap()
    .unwrap();
    let abundance = [2_usize, 3, 5];
    for iteration in 0..3 {
        observer
            .observe(TrajectoryObservation {
                iteration,
                physical_time: None,
                abundance: AbundanceView::Counts(&abundance),
                detector_observable: None,
                equilibrium_evidence: EquilibriumEvidence::Unavailable,
            })
            .unwrap();
    }
    let before = ALLOCATIONS.load(Ordering::Relaxed);
    let signal = observer
        .observe(TrajectoryObservation {
            iteration: 3,
            physical_time: None,
            abundance: AbundanceView::Counts(&abundance),
            detector_observable: None,
            equilibrium_evidence: EquilibriumEvidence::Unavailable,
        })
        .unwrap();
    let after = ALLOCATIONS.load(Ordering::Relaxed);
    assert!(signal.is_none());
    assert_eq!(after, before);
}
