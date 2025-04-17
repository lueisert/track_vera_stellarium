# Track Vera C. Rubin Observatory in Stellarium

Script for using the *Mosaic Camera* Stellarium plugin with Vera C. Rubin Observatory.


Mosaic Camera Plugin
====================

To get the plugin, you will need a recent version of [Stellarium](https://stellarium.org/).  Until Version 25.2 is released in Summer 2025, that means downloading a [weekly snapshot](https://github.com/Stellarium/stellarium-data/releases/tag/weekly-snapshot).

To enable the Mosaic Camera plugin, first load the Configurations window (F2 by default), and locate the Plugins tab.  Scroll down to Mosaic Camera and check the "Load at startup" box.  Restart Stellarium.

Once you've restarted, you can navigate back to the Mosaic Camera plugin and click the "configure" button, which will bring up a dialog for different cameras and their current RA, dec, rotation, and visibility parameters.  (You can also right-click on the Mosaic Camera icon that appeared in the toolbar).  From the dialog you can also recenter the Stellarium view on the current camera, or set the camera position to the currently selected object or current view.

It's possible to associate many of these functions with keyboard shortcuts too.  (Press F7 to bring up the Keyboard Shortcuts editor.)


tvs.py
======
The Mosaic Camera plugin is also controllable over an http connection, which is where the tvs.py script comes in.  With this script you can:
- slew the camera mosaic overlay to a given RA/DEC/RotSkyPos or Az/Alt/RotTelPos
- slew to a target by name
- slew to a specific past visit (requires access to the Rubin Consolidated Database)
- replay a past night's observations (requires access to the Rubin Consolidated Database)
- follow the Rubin (and LATISS) cameras in real time (accesses the Rubin Engineering Facilities Database)
- display detailed status of the current camera overlay

Note that for these functions to work optimally, you should set your Stellarium location to Cerro Pachón (Use F6).

Requirements are `click` and `astropy` for everything.  You'll need `astroquery` to slew to a target by name.  If you have the LSST Science Pipelines `summit_utils` package available, that'll be used to access the Rubin Consolidated Database (but I've added a rudimentary fallback if this is unavailable.)  Similarly, if `lsst_efd_client` is available, that'll be used to access the Rubin Engineering Facilities Database, but a rudimentary fallback also exists.

There's no installation.  The script is self-contained.  You can either call it with python a la `python tvs.py ...` or add your own [shebang](https://en.wikipedia.org/wiki/Shebang_(Unix)) to the beginning and just use `tvs.py ...`.


Commands
========
tvs.py works with subcommands, kinda like `git`.

Running

`python tvs.py`

Will display the overall help message including listing available subcommands.


python tvs.py slew
------------------
The slew command will reposition and set visible a camera mosaic overlay.  For example:

`python tvs.py slew 12:34:56.7 -12:34:56.7 45`

Slew the LSSTCam mosaic overlay to RA 12h 34m 56.7s, Dec -12d 34m 56.7s and sets the sky rotation angle to 45 degrees (angle of the +Y axis in data visualization coordinate system (DVCS) measured from North through East).

`python tvs.py slew 12:34:56.7 -12:34:56.7 45 --horizon`

Slew the LSSTCam mosaic overlay to azimuth 12 _degrees_ 34m 56.7s, altitude -12d 34m 56.7s and sets the hardware rotation angle to 45 degrees (angle of -Y axis is the camera coordinate system (CCS) measured from +alt toward +az).

Note you can use other strings that `astropy.coordinates.Angle` can interpret, e.g.,

`python tvs.py slew 12h34m56.7 -12d34m56.7 45deg`

`python tvs.py slew 12.5h -12.5d "45 deg"`

To slew a camera other than LSSTCam, use the `--camera` argument, e.g.,

`python tvs.py slew "185 deg" "12 deg" "0.5 rad" --camera LATISS`

By default, the slew and some other commands will change the view within Stellarium to be centered around the newly positioned overlay.  To keep the view the same, use the `--no-follow` flag:

`python tvs.py slew 12:34:56.7 -12:34:56.7 45 --no-follow`

To also set the time inside of Stellarium, you can use

`python tvs.py slew 0 0 0 --time 2024-04-08T00:00:00`

`python tvs.py slew 0 0 0 --time 60408.0 --timeformat mjd`


python tvs.py target
--------------------
Pretty much the same as slew, but instead of providing spherical coords, provide a name to lookup via `astroquery`.  You can still set the rotator parameter as in `slew`.  Some examples:

Keep current RotSkyPos

`python tvs.py target M20`

Set RotSkyPos to 45 degrees

`python tvs.py target M20 45`

Set RotTelPos to 45 degrees

`python tvs.py target Trifid_Nebula 45 --horizon`

And so on:

`python tvs.py target M20 20deg`

`python tvs.py target M20 "20 deg"`

`python tvs.py target "Trifid Nebula" 45deg --camera LATISS`

`python tvs.py target M20 45 --no-follow`

`python tvs.py target M20 --time 2024-04-08T00:00:00`

`python tvs.py target M20 --time 60408.0 --timeformat mjd`


python tvs.py visit
-------------------
The `visit` subcommand will look up a past visit in the Rubin Consolidated Database (ConsDB).  It will attempt to infer the correct camera based on the visit number, but you can also explicitly set one of ComCam, LSSTCam, or LATISS.

To access ConsDB, you must have a Rubin Science Platform security token.  (See https://nb.lsst.io/environment/tokens.html#using-a-token-outside-the-science-platform).  tvs.py will look for a token at `~/.lsst/rsp_token` by default, but you can specify a path explicitly with `--rsp-token`.

Some examples:

Infers ComCam from date 2024-12-11

`python tvs.py visit 2024121100608`

Infers LATISS from date 2025-04-03

`python tvs.py visit 2025040300825`

Specify camera explicitly

`python tvs.py visit 2024121100608 --camera ComCam`

`python tvs.py visit 2025040300825 --camera LATISS`

`python tvs.py visit 2025040300825 --camera LATISS --no-follow`


python tvs.py replay
--------------------
The `replay` subcommand takes a date (and optional camera), and uses the Stellarium internal clock to overlay the appropriate mosaic.  Like `visit`, this subcommand requires ConsDB access.  Some examples:

Replay night of LSSTCam first photon:

`python tvs.py replay 20250415`

Replay last two nights of ComCam:

`python tvs.py replay 20241210 20241211`


python tvs.py follow
--------------------
The `follow` subcommand polls the Rubin Engineering Facilities Database (EFD) for the most recent pointing and rotation of either LSSTCam or LATISS to use to update the stellarium mosaic overlays.  The EFD is currently world public (though that may change) so no special credentials are required.  Example:

`python tvs.py follow`

`python tvs.py follow LATISS`

If you want to follow both LSSTCam and LATISS, you can run in separate terminals both of the above (though you'll probably want to use `--no-follow` on at least one of them so the Stellarium view isn't constantly flipping back and forth between the two telescopes.)


python tvs.py status
--------------------
This command simply prints to stdout the current position (in both equatorial and horizon coordinates) of the specified camera overlay.  Example:

`python tvs.py status`

`python tvs.py status LATISS`
