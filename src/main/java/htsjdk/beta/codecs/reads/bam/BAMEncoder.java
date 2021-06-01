package htsjdk.beta.codecs.reads.bam;

import htsjdk.beta.plugin.bundle.Bundle;
import htsjdk.beta.plugin.bundle.BundleResourceType;
import htsjdk.beta.plugin.reads.ReadsEncoderOptions;
import htsjdk.beta.plugin.reads.ReadsFormat;
import htsjdk.beta.plugin.reads.ReadsEncoder;

// TODO: handle presorted

/**
 * Base class for BAM encoders.
 */
public abstract class BAMEncoder implements ReadsEncoder {
    protected final Bundle outputBundle;
    protected final ReadsEncoderOptions readsEncoderOptions;
    final private String displayName;

    public BAMEncoder(final Bundle outputBundle, final ReadsEncoderOptions readsEncoderOptions) {
        this.outputBundle = outputBundle;
        this.readsEncoderOptions = readsEncoderOptions;
        this.displayName = outputBundle.get(BundleResourceType.READS).get().getDisplayName();
    }

    @Override
    final public ReadsFormat getFormat() { return ReadsFormat.BAM; }

    @Override
    final public String getDisplayName() { return displayName; }

}