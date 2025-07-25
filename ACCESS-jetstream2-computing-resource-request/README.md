# Application submission info

---
> Submitted on 15-Apr-2025; approved 17-Apr-2025 (BIO250131)
---

Submitted through ACCESS: https://access-ci.org/

Submitted an "Explore ACCESS", different types listed here: https://allocations.access-ci.org/project-types

> In the years prior, this has been the way to go. In 2025, it became a huge ordeal to get teh full 400,000 credits up front. So in the future, submit a "Discover ACCESS" request instead.

After logged in, submission process started from this page: https://allocations.access-ci.org/opportunities

Then selecting to "Request New Project", and choosing "Request a Explore ACCESS project".

## Required information I added to the form

**Title:** STAMPS 2025 at Marine Biologial Laboratory in Woods Hole, MA, USA

**Public overview:**  

> General overview
> The STAMPS course (Strategies and Techniques for Analyzing Microbial Population Structures) has been a yearly event at the Marine Biological Laboratory (MBL) in Woods Hole, MA, USA for over a decade, helping nearly 1,000 learners establish a foundation in bioinformatics over the years. Since 2019 we have been fortunate enough to have our computational infrastructure provided through XSEDE, then ACCESS, and Jetstream (1 and 2), and we hope to continue that this year. The course takes place this year 14-July to 24-July and will involve a lot of command-line, R, amplicon analysis, metagenomics, statistical training, and other concepts and computational activities (last year's schedule can be seen here: https://github.com/mblstamps/stamps2024/wiki#schedule).
> 
> How we plan to use ACCESS resources
> We will be running an about 11-day workshop (14-July to 24-July). I have experience in the past doing this through Indiana Jetstream2, managing through Exosphere, so that would be great if possible. We are requesting the ability to run 60 (accounting for faculty and participants) m3.large instances concurrently for about 15 days straight. This is longer than we plan on them being active, but would like to leave room for building the image and testing/setting things up ahead of time, and any unexpected situations. Using the Jetstream2 estimator (https://docs.jetstream-cloud.org/alloc/estimator/), with 60 m3.large instances for 15 days, this comes out to a request for 345,600 SUs/ACCESS Credits. We think this will work under the 400,000 max for "Explore ACCESS" allocations. Since we will be using these all at once over a short time, we do request the full amount to be released up front if possible. Please let me know if I should be estimating a different way or if any other information is needed.
> 
> Thanks for your consideration and any help!

**Keywords:** microbial ecology, bioinformatics

**How do you plan to use this project?** Classroom / Training

**Opportunity questions:**
- Would you like us to provide an advisory review? No
- Is the planned work associated with any of the following? No
- How did you hear about ACCESS? Word of mouth

**Fields of Science:** “Other Biological Sciences”

**Supporting Grants?** No

**Documents:** Needed to attach my CV here.

**Available resources:** checked box for ACCESS Credits

---

## After submission

### Transferring credits from ACCESS to JetStream2

Once approved, and logged in, needed this page (https://allocations.access-ci.org/requests) in order to transfer ACCESS credits to a specific resource. For the appropriate allocation/Project, selected "Credits + Resources", then the text box that initially says "Add a resource to your exchange...", then selected "Indiana Jetstream2 CPU", then entered 196,000 (since I did this before the other 200,000 were available). Then for "Indiana Jetstream2 Storage" added the remaining 4,000, giving 4 TB of shared storage to use.
After the other 200,000 are made available (see next section), repeat the process without adding anymore to Storage.

### Getting the remaining 200,000 credits

The “Explore ACCESS” project at this time comes with a total of 400,000 credits, with 200,000 released at first. After getting approved, I went to the https://allocations.access-ci.org/requests page, selected the STAMPS 2025 project, clicked on “Credits + Resources”, then “REQUEST MORE CREDITS”, then “REQUEST MORE CREDITS” again, then “REQUEST A SUPPLEMENT”. Then filled out these:

**Reason**
> We have not used any credits for this allocation yet. This allocation is for a course that will be using all the credits (nearly 400,000) rapidly over a 10-11 day period. In the past few years, we have been granted the full 400,000 right up front, which is extremely helpful as it prevents jeopardizing the 60+ person course while allowing credits to run low prior to requesting the supplement. We have learned this year that this is not ideal, and in the future we will request a "Discover ACCESS" allocation from the start to avoid added problems like this. But for this instance, could we please have the remaining 200,000 credits now?
> 
> Sorry for the complications, and thanks!
> -Mike


- Available Resources:
  -	Checked the box for ACCESS Credits
- Document Type:
  -	A Progress Report is required. This was just a PDF stating the same as the “Reason” above.

Then submitted the form.

### Requesting quota-limit increases so we can run up to 60 instances concurrently, and requesting Manila (for shared volume)

The allocation comes with a limit on the number of concurrent instances that can be run. I submitted a request when logged into Jetstream2 here: https://jetstream2.exosphere.app/exosphere/getsupport

This is the text I submitted (after selecting the radio dial for “An Allocation”):

> Hi there :)
> 
> We plan to use this allocation (BIO250131) with 60 concurrent m3.large instances for a bioinformatics course we are running. But the starting limits prohibit that.
> 
> Could you please help with increasing the required quotas so that we will be able to run up to 60 m3.large instances concurrently on this allocation, including cores, ram, volume, ports, available IP addresses, and whatever other magic you folks take care of?
> 
> We would also like to be able to use the 4 TB we specified for Jetstream2 Storage as a shared volume attached to all instances. We've done this with Manila in the past, managed through Horizon. Could you please enable that for us/help with whatever is needed for me to be able to set that up? Ultimately, I plan to create the base instance with the share mounted, then imaging that to use for creating the instances for the course.
> 
> Thank you for any help!  
> -Mike

And actually, that site just helps make the email you can then copy/paste to send them via regular email. 

## Instance creation and setup info
Instance creation and setup info for the 2024 year is here: https://hackmd.io/kIeWoqZYSJCltfkeLjfDtA?view
