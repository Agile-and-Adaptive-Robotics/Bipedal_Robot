// Question slides only. Original app screenshots remain unchanged on disk.
// Each crop omits the lower selected-record pane instead of altering its text.
export async function apply(h) {
  const {begin, im, tx, old, GREEN, INK} = h;
  const source = 'Dissertation Ch. 5 · Sensory Afferent Database (SADb)';

  begin(117, 'SADb connects a selected paper to its references and citing papers',
    'SADb is the Sensory Afferent Database. These are enlarged, authentic portions of the selected-record and connected-paper panes. The selected Grillner (1985) record retains its bibliographic identity, DOI, species tags, and graph-navigation actions. The adjacent pane lists connected papers. The crop stops before the record curation text. Citation relationships are literature metadata, not biological synapses. Source: submitted Chapter 5 and the local SADb app. Example paper: Grillner, S. (1985), Neurobiological Bases of Rhythmic Motor Acts in Vertebrates, Science, DOI 10.1126/science.3975635.',
    source);
  im(117, old('s117_img0.png'), 110, 126, 300, 315,
    {left: 0, top: 0, right: 0, bottom: .46});
  im(117, old('s117_img1.png'), 550, 126, 300, 315,
    {left: 0, top: 0, right: 0, bottom: .46});
  tx(117, 'Selected record', 70, 445, 370, 32, 25, true, GREEN, 'center');
  tx(117, 'Connected papers', 520, 445, 370, 32, 25, true, GREEN, 'center');

  begin(118, 'SADb combines searchable records with an interactive citation map',
    'This authentic app screenshot shows the Table, Pivot, Bubble map, and Search tabs, the map controls, the selected record, and connected-paper navigation. The lower curation pane is outside the crop. Local Table, Pivot, and Bubble map views work offline. Search (online), DOI destinations, and external papers require an internet connection. The full collection contains 1,393 literature records; 943 is the Papers category displayed by this app snapshot. It is not a deduplicated count of the full collection. Source: submitted Chapter 5 and the local SADb app.',
    source);
  im(118, old('s123_img0.png'), 30, 115, 900, 328,
    {left: 0, top: 0, right: 0, bottom: .27});
  tx(118, 'Select a paper to explore its citation neighborhood',
    45, 447, 870, 32, 25, true, GREEN, 'center');

  // Source 119 already shows the unaltered full-table view and 1,393 / 943 distinction.

  begin(120, 'A curated finding can suggest a testable robot hypothesis',
    'The left statement is a concise paraphrase of the Grillner (1985) record, not a verbatim quotation. The right question is a proposed robot or simulation experiment recorded separately from the source finding. It has not been performed for this example. A useful test would compare shared and separate subnetworks while holding the mechanical model and input protocol fixed, then measure motor-pattern timing and stability. The database preserves the distinction between what the literature reports and what a researcher proposes to test. Source: submitted Chapter 5; Grillner (1985), Neurobiological Bases of Rhythmic Motor Acts in Vertebrates, Science, DOI 10.1126/science.3975635.',
    'SADb curation example · Grillner (1985) · DOI: 10.1126/science.3975635');
  tx(120, 'Literature finding', 50, 133, 390, 38, 27, true, GREEN);
  tx(120, 'Spinal networks generate rhythmic motor patterns.',
    50, 202, 390, 170, 30, false, INK);
  tx(120, 'Proposed robot test', 520, 133, 390, 38, 27, true, GREEN);
  tx(120, 'Can a shared neural subnetwork preserve timing across two simulated motor patterns?',
    520, 202, 390, 215, 30, false, INK);

  begin(121, 'Citation links help researchers follow related work',
    'The authentic Bubble map view focuses the literature graph around Grillner (1985). The selected-record pane is outside this crop. Colored direct links distinguish references from citing papers; gray links extend the displayed citation neighborhood. The Connections pane supports navigation to the related records. The neuron-like display is a visual metaphor for a citation network, not an inferred neural circuit. Source: submitted Chapter 5 and the local SADb app.',
    source);
  im(121, old('s121_img0.png'), 45, 111, 870, 337,
    {left: .216, top: .173, right: 0, bottom: 0});
  tx(121, 'Red: references     Green: citing papers',
    50, 449, 860, 31, 24, true, GREEN, 'center');

  // Source 122 already shows the original pivot view, including its counting rule.

  // Source 123 owns the authentic embedded offline walkthrough, not just a poster.
  // Preserve its full body, media relationship, playback timing, hyperlink, and
  // existing notes. It must not be marked changed or transplanted from the imported
  // artifact representation, which may flatten the video to its poster image.
  // Its original OOXML part in v19 is ppt/slides/slide28.xml, not slide123.xml.
}
